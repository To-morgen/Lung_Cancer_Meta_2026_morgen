#!/usr/bin/env Rscript
# ============================================================================
# 02_cluster_qc.R — Cluster-level QC & composition analysis
#
# Generates:
#   1. Clustree (resolution tree) — if clustree package available
#   2. Cluster composition by sample / group (bar charts)
#   3. Per-cluster QC metrics (nFeature, nCount, MT%, cell cycle)
#   4. Cluster size distribution
#   5. Summary tables for Notion/review
#
# Key decision point: review plots → choose final resolution →
#   update scrna_qc_params.yaml: clustering.default_resolution
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(data.table)
})

source(here("workflow", "scrna", "functions", "io_scrna.R"))
source(here("workflow", "scrna", "functions", "qc_utils.R"))
source(here("scripts", "utils", "utils_plotting.R"))

# ---- Config ----
QC_PARAMS  <- load_qc_params()
clust_dirs <- scrna_cluster_dirs()

RESOLUTIONS <- QC_PARAMS$clustering$resolutions %||% c(0.3, 0.5, 0.8, 1.0, 1.2)
DEFAULT_RES <- QC_PARAMS$clustering$default_resolution %||% 0.8

cat("\n╔══════════════════════════════════════════════════╗\n")
cat("║   Phase 5 Step 02: Cluster QC & Composition      ║\n")
cat(sprintf("║   Default resolution: %.1f                        ║\n", DEFAULT_RES))
cat("╚══════════════════════════════════════════════════╝\n\n")

# ---- Load ----
sobj_file <- file.path(clust_dirs$objects, "seurat_clustered.rds")
if (!file.exists(sobj_file)) stop("Clustered object not found: ", sobj_file)

log_msg(sprintf("Loading: %s", basename(sobj_file)))
sobj <- readRDS(sobj_file)
md <- sobj@meta.data

default_col <- sprintf("clusters_res%.1f", DEFAULT_RES)
if (!default_col %in% colnames(md)) {
  stop(sprintf("Column %s not found. Available: %s",
               default_col,
               paste(grep("clusters_res", colnames(md), value = TRUE), collapse = ", ")))
}

log_msg(sprintf("Using resolution %.1f → %d clusters",
                DEFAULT_RES, length(unique(md[[default_col]]))))

# ============================================================================
# Plot 1: Clustree (resolution stability)
# ============================================================================
log_msg("Plot 1: Clustree...")

clustree_available <- requireNamespace("clustree", quietly = TRUE)
if (clustree_available) {
  library(clustree)

  # clustree needs columns named "res.X"
  # Seurat already creates SCT_snn_res.X columns
  res_cols <- grep("^SCT_snn_res\\.", colnames(md), value = TRUE)

  if (length(res_cols) >= 2) {
    pdf(file.path(clust_dirs$plots, "04_clustree.pdf"), width = 12, height = 10)
    print(clustree(sobj, prefix = "SCT_snn_res.") +
            labs(title = "Clustree: Resolution Stability") +
            theme(plot.title = element_text(face = "bold", size = 14)))
    dev.off()
    png(file.path(clust_dirs$plots, "04_clustree.png"), width = 1200, height = 1000, res = 150)
    print(clustree(sobj, prefix = "SCT_snn_res.") +
            labs(title = "Clustree: Resolution Stability"))
    dev.off()
    log_msg("  ✅ 04_clustree saved")
  }
} else {
  log_msg("  ⚠️ clustree not installed — skipping. Install with: install.packages('clustree')")
}

# ============================================================================
# Plot 2: Cluster composition by sample & group
# ============================================================================
log_msg("Plot 2: Cluster composition...")

md$cluster <- md[[default_col]]

# Cell counts per cluster per sample
comp_sample <- md %>%
  count(cluster, sample_id) %>%
  group_by(sample_id) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p_comp1 <- ggplot(comp_sample, aes(x = cluster, y = pct, fill = sample_id)) +
  geom_col(position = "dodge", alpha = 0.8) +
  labs(title = "Cluster Composition by Sample (%)", x = "Cluster", y = "% of sample", fill = "Sample") +
  theme_project() +
  theme(axis.text.x = element_text(angle = 0))

# Cell proportions per cluster per group
comp_group <- md %>%
  count(cluster, group) %>%
  group_by(group) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p_comp2 <- ggplot(comp_group, aes(x = cluster, y = pct, fill = group)) +
  geom_col(position = "dodge", alpha = 0.8) +
  labs(title = "Cluster Composition by Group (%)", x = "Cluster", y = "% of group", fill = "Group") +
  theme_project()

# Stacked bar: what samples contribute to each cluster
comp_cluster <- md %>%
  count(cluster, sample_id) %>%
  group_by(cluster) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

p_comp3 <- ggplot(comp_cluster, aes(x = cluster, y = pct, fill = sample_id)) +
  geom_col(alpha = 0.8) +
  labs(title = "Sample Contribution per Cluster (stacked %)",
       x = "Cluster", y = "% of cluster", fill = "Sample") +
  theme_project()

pdf(file.path(clust_dirs$plots, "05_cluster_composition.pdf"), width = 16, height = 16)
print(p_comp1 / p_comp2 / p_comp3)
dev.off()
png(file.path(clust_dirs$plots, "05_cluster_composition.png"), width = 1600, height = 1600, res = 150)
print(p_comp1 / p_comp2 / p_comp3)
dev.off()
log_msg("  ✅ 05_cluster_composition saved")

# Save composition table
fwrite(comp_sample, file.path(clust_dirs$reports, "cluster_composition_by_sample.csv"))
fwrite(comp_group, file.path(clust_dirs$reports, "cluster_composition_by_group.csv"))

# ============================================================================
# Plot 3: Per-cluster QC metrics
# ============================================================================
log_msg("Plot 3: Per-cluster QC metrics...")

qc_metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
if ("percent.ribo" %in% colnames(md)) qc_metrics <- c(qc_metrics, "percent.ribo")
if ("S.Score" %in% colnames(md)) qc_metrics <- c(qc_metrics, "S.Score", "G2M.Score")

qc_plots <- list()
for (metric in qc_metrics) {
  qc_plots[[metric]] <- VlnPlot(sobj, features = metric, group.by = default_col,
                                  pt.size = 0, log = (metric %in% c("nCount_RNA"))) +
    labs(title = metric) +
    theme_project() + NoLegend() +
    theme(axis.text.x = element_text(size = 7))
}

pdf(file.path(clust_dirs$plots, "06_cluster_qc_metrics.pdf"), width = 18, height = 12)
print(wrap_plots(qc_plots, ncol = 3) +
        plot_annotation(title = sprintf("QC Metrics per Cluster (res=%.1f)", DEFAULT_RES),
                        theme = theme(plot.title = element_text(face = "bold", size = 14))))
dev.off()
png(file.path(clust_dirs$plots, "06_cluster_qc_metrics.png"), width = 1800, height = 1200, res = 150)
print(wrap_plots(qc_plots, ncol = 3) +
        plot_annotation(title = sprintf("QC Metrics per Cluster (res=%.1f)", DEFAULT_RES)))
dev.off()
log_msg("  ✅ 06_cluster_qc_metrics saved")

# ============================================================================
# Plot 4: Cluster sizes
# ============================================================================
log_msg("Plot 4: Cluster sizes...")

size_df <- md %>%
  count(cluster) %>%
  arrange(desc(n)) %>%
  mutate(pct = round(n / sum(n) * 100, 1),
         cluster = factor(cluster, levels = cluster))

p_size <- ggplot(size_df, aes(x = cluster, y = n, fill = cluster)) +
  geom_col(alpha = 0.8, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%s\n(%.1f%%)", scales::comma(n), pct)),
            vjust = -0.3, size = 2.5) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  labs(title = sprintf("Cluster Sizes (res=%.1f, %d clusters, %s cells total)",
                       DEFAULT_RES, nrow(size_df), scales::comma(sum(size_df$n))),
       x = "Cluster", y = "Cells") +
  theme_project()

pdf(file.path(clust_dirs$plots, "07_cluster_sizes.pdf"), width = 14, height = 6)
print(p_size)
dev.off()
png(file.path(clust_dirs$plots, "07_cluster_sizes.png"), width = 1400, height = 600, res = 150)
print(p_size)
dev.off()
log_msg("  ✅ 07_cluster_sizes saved")

# Save size table
fwrite(size_df, file.path(clust_dirs$reports, "cluster_sizes.csv"))

# ============================================================================
# Summary table for all resolutions
# ============================================================================
log_msg("Generating resolution comparison summary...")

res_compare <- data.frame()
for (res in RESOLUTIONS) {
  col <- sprintf("clusters_res%.1f", res)
  if (col %in% colnames(md)) {
    clusters <- md[[col]]
    n_cl <- length(unique(clusters))
    sizes <- table(clusters)
    res_compare <- rbind(res_compare, data.frame(
      resolution     = res,
      n_clusters     = n_cl,
      min_size       = min(sizes),
      max_size       = max(sizes),
      median_size    = median(sizes),
      smallest_pct   = round(min(sizes) / nrow(md) * 100, 2),
      is_default     = ifelse(res == DEFAULT_RES, "<<<", ""),
      stringsAsFactors = FALSE
    ))
  }
}

fwrite(res_compare, file.path(clust_dirs$reports, "resolution_comparison.csv"))

cat("\n=== Resolution Comparison ===\n")
print(res_compare)

cat(sprintf("\n💡 Review plots in: %s", clust_dirs$plots))
cat(sprintf("\n💡 Current default: res=%.1f", DEFAULT_RES))
cat("\n💡 To change: edit configs/params/scrna_qc_params.yaml → clustering.default_resolution")
cat("\n\n========== Cluster QC complete ==========\n")
