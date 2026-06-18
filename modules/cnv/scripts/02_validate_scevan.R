#!/usr/bin/env Rscript
# ============================================================================
# 02_validate_scevan.R — Validate SCEVAN recovery results
#
# Checks:
#   1. Tumor/non-tumor label distribution by cluster, sample, group
#   2. CNA signal heatmap (from CNA_mtx_relat)
#   3. Tumor fraction per group (FL vs A1 vs mc)
#   4. Known biology sanity checks (ref clusters should be ~0% tumor)
#   5. UMAP overlay of tumor labels
#
# Input:  seurat_with_scevan.rds (from 03_recover)
#         SCEVAN .RData artifacts
# Output: validation plots + CSV reports
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(yaml)
})

# ---- Paths ----
PROJ_ROOT <- Sys.getenv("LUNGMETA_ROOT", unset = normalizePath(here::here("..", "..")))
module_root <- here::here()

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

path_from_root <- function(path) {
  if (grepl("^/", path)) path else file.path(PROJ_ROOT, path)
}

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else if (nzchar(Sys.getenv("CNV_CONFIG", unset = ""))) {
  Sys.getenv("CNV_CONFIG")
} else {
  file.path(module_root, "configs", "cnv_params.yaml")
}
if (!file.exists(cfg_path)) stop("Config not found: ", cfg_path)
cfg <- yaml::read_yaml(cfg_path)

SCEVAN_DIR <- if (!is.null(cfg$output$scevan)) {
  path_from_root(cfg$output$scevan)
} else {
  file.path(path_from_root(cfg$output$main_results), "scevan")
}
REPORT_DIR <- if (!is.null(cfg$output$reports)) {
  path_from_root(cfg$output$reports)
} else {
  file.path(path_from_root(cfg$output$main_results), "reports")
}
PLOT_DIR <- if (!is.null(cfg$output$plots)) {
  path_from_root(cfg$output$plots)
} else {
  file.path(path_from_root(cfg$output$main_results), "plots")
}
SCEVAN_OUT <- Sys.getenv("CNV_LOCAL_SCEVAN_DIR", unset = "")
if (SCEVAN_OUT == "") {
  SCEVAN_OUT <- if (!is.null(cfg$output$module_scevan)) {
    path_from_root(cfg$output$module_scevan)
  } else {
    file.path(module_root, cfg$output$module_results %||% "results", "scevan")
  }
}

for (d in c(REPORT_DIR, PLOT_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Config ----
ref_clusters <- as.character(cfg$scevan$reference_clusters %||% character())

cat("\n")
cat("╔══════════════════════════════════════════════════╗\n")
cat("║   02: SCEVAN Validation                          ║\n")
cat(sprintf("║   Time: %s              ║\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("╚══════════════════════════════════════════════════╝\n\n")

# ============================================================================
# 1. Load recovered Seurat
# ============================================================================
cat("Loading recovered Seurat...\n")
sobj_path <- file.path(SCEVAN_DIR, "seurat_with_scevan.rds")
if (!file.exists(sobj_path)) stop("Recovered Seurat not found: ", sobj_path)
sobj <- readRDS(sobj_path)

cat(sprintf("  Cells: %d\n", ncol(sobj)))
cat(sprintf("  Metadata cols: %s\n", paste(colnames(sobj@meta.data), collapse = ", ")))


# Check scevan_label exists
if (!"scevan_label" %in% colnames(sobj@meta.data)) {
  stop("scevan_label not in metadata — run 03_recover first")
}

md <- sobj@meta.data

cat("\n=== Overall Label Distribution ===\n")
print(table(md$scevan_label))
cat("\n")

# ============================================================================
# 2. Tumor fraction by cluster
# ============================================================================
cat("=== Tumor Fraction by Cluster ===\n")

cluster_stats <- md %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_total    = n(),
    n_tumor    = sum(scevan_label == "tumor", na.rm = TRUE),
    n_nontumor = sum(scevan_label == "non_tumor", na.rm = TRUE),
    pct_tumor  = round(n_tumor / n() * 100, 1),
    .groups    = "drop"
  ) %>%
  arrange(desc(pct_tumor))

cluster_stats$is_ref <- as.character(cluster_stats$seurat_clusters) %in% ref_clusters

print(as.data.frame(cluster_stats), row.names = FALSE)
fwrite(cluster_stats, file.path(REPORT_DIR, "scevan_tumor_fraction_by_cluster.csv"))

# Sanity: ref clusters should be ~0% tumor
ref_tumor <- cluster_stats %>% filter(is_ref == TRUE)
if (nrow(ref_tumor) > 0) {
  cat(sprintf("\n  Ref cluster tumor %%: min=%.1f%%, max=%.1f%%, mean=%.1f%%\n",
              min(ref_tumor$pct_tumor), max(ref_tumor$pct_tumor), mean(ref_tumor$pct_tumor)))

  if (max(ref_tumor$pct_tumor) > 5) {
    cat("  ⚠️  WARNING: Some ref clusters have >5% tumor — check ref selection\n")
  } else {
    cat("  ✅ Ref clusters are clean (<5% tumor)\n")
  }
} else {
  cat("\n  ⚠️  No configured ref clusters found in SCEVAN output\n")
}

# ============================================================================
# 3. Tumor fraction by sample & group
# ============================================================================
cat("\n=== Tumor Fraction by Sample ===\n")

sample_stats <- md %>%
  group_by(sample_id, group) %>%
  summarise(
    n_total   = n(),
    n_tumor   = sum(scevan_label == "tumor", na.rm = TRUE),
    pct_tumor = round(n_tumor / n() * 100, 1),
    .groups   = "drop"
  ) %>%
  arrange(group, sample_id)

print(as.data.frame(sample_stats), row.names = FALSE)
fwrite(sample_stats, file.path(REPORT_DIR, "scevan_tumor_fraction_by_sample.csv"))

cat("\n=== Tumor Fraction by Group ===\n")

group_stats <- md %>%
  group_by(group) %>%
  summarise(
    n_total   = n(),
    n_tumor   = sum(scevan_label == "tumor", na.rm = TRUE),
    pct_tumor = round(n_tumor / n() * 100, 1),
    .groups   = "drop"
  )

print(as.data.frame(group_stats), row.names = FALSE)
fwrite(group_stats, file.path(REPORT_DIR, "scevan_tumor_fraction_by_group.csv"))

# ============================================================================
# 4. Identify tumor-enriched clusters (candidate tumor clusters)
# ============================================================================
cat("\n=== Candidate Tumor Clusters (>50% tumor) ===\n")
tumor_clusters <- cluster_stats %>%
  filter(pct_tumor > 50) %>%
  pull(seurat_clusters)

cat(sprintf("  Clusters: %s\n", paste(tumor_clusters, collapse = ", ")))

mixed_clusters <- cluster_stats %>%
  filter(pct_tumor > 10 & pct_tumor <= 50) %>%
  pull(seurat_clusters)

if (length(mixed_clusters) > 0) {
  cat(sprintf("  Mixed clusters (10-50%%): %s\n", paste(mixed_clusters, collapse = ", ")))
}

# ============================================================================
# 5. Plots
# ============================================================================
cat("\nGenerating validation plots...\n")

# ---- Plot 1: UMAP by scevan_label ----
tumor_cols <- c("tumor" = "#E41A1C", "non_tumor" = "#377EB8", "ref_normal" = "#999999")

p1 <- DimPlot(sobj, group.by = "scevan_label", cols = tumor_cols, pt.size = 0.1) +
  labs(title = "SCEVAN Recovery: Tumor vs Non-tumor") +
  theme(plot.title = element_text(face = "bold"))

# ---- Plot 2: Tumor fraction bar per cluster ----
p2_data <- cluster_stats %>%
  mutate(cluster_label = paste0("C", seurat_clusters,
                                 ifelse(is_ref, " (ref)", "")))

p2 <- ggplot(p2_data, aes(x = reorder(cluster_label, -pct_tumor), y = pct_tumor)) +
  geom_col(aes(fill = pct_tumor > 50), alpha = 0.8) +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
                    labels = c("TRUE" = ">50% tumor", "FALSE" = "<=50% tumor")) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red", alpha = 0.5) +
  labs(title = "Tumor Fraction per Cluster", x = "", y = "% Tumor", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# ---- Plot 3: Tumor fraction by sample ----
p3 <- ggplot(sample_stats, aes(x = sample_id, y = pct_tumor, fill = group)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", pct_tumor)), vjust = -0.3, size = 3) +
  labs(title = "Tumor Fraction per Sample", x = "", y = "% Tumor") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- Plot 4: UMAP split by group, colored by scevan_label ----
p4 <- DimPlot(sobj, group.by = "scevan_label", split.by = "group",
              cols = tumor_cols, pt.size = 0.05, ncol = 3) +
  labs(title = "SCEVAN Labels by Group")

# ---- Save composite ----
pdf(file.path(PLOT_DIR, "03_scevan_validation.pdf"), width = 18, height = 18)
print((p1 | p2) / (p3 | plot_spacer()) / p4 + plot_layout(heights = c(1, 0.8, 1)))
dev.off()

png(file.path(PLOT_DIR, "03_scevan_validation.png"), width = 1800, height = 1800, res = 150)
print((p1 | p2) / (p3 | plot_spacer()) / p4 + plot_layout(heights = c(1, 0.8, 1)))
dev.off()

cat("  ✅ 03_scevan_validation\n")

# ============================================================================
# 6. CNA Heatmap (if CNA_mtx_relat available)
# ============================================================================
cat("\nCNA Heatmap...\n")

# Find the .RData file
rdata_files <- list.files(SCEVAN_OUT, pattern = "\\.RData$", full.names = TRUE, recursive = TRUE)

if (length(rdata_files) > 0) {
  rdata_file <- rdata_files[1]
  cat(sprintf("  Loading: %s\n", basename(rdata_file)))

  env <- new.env()
  load(rdata_file, envir = env)

  if ("CNA_mtx_relat" %in% ls(env)) {
    cna_mtx <- env$CNA_mtx_relat
    cat(sprintf("  CNA matrix: %d genes × %d cells\n", nrow(cna_mtx), ncol(cna_mtx)))

    # Match barcodes to Seurat
    common_bc <- intersect(colnames(cna_mtx), colnames(sobj))
    cat(sprintf("  Matched to Seurat: %d cells\n", length(common_bc)))

    if (length(common_bc) > 500) {
      # Subsample for heatmap (max 3000 cells, balanced)
      set.seed(42)
      md_sub <- md[common_bc, ]

      max_per_label <- 1500
      tumor_bc <- rownames(md_sub[md_sub$scevan_label == "tumor", ])
      nontumor_bc <- rownames(md_sub[md_sub$scevan_label != "tumor", ])

      if (length(tumor_bc) > max_per_label) tumor_bc <- sample(tumor_bc, max_per_label)
      if (length(nontumor_bc) > max_per_label) nontumor_bc <- sample(nontumor_bc, max_per_label)
      sel_bc <- c(tumor_bc, nontumor_bc)

      cna_sub <- as.matrix(cna_mtx[, sel_bc])

      # Clip extreme values
      cna_sub[cna_sub > 2] <- 2
      cna_sub[cna_sub < -2] <- -2

      # Annotation
      anno_df <- data.frame(
        Label   = md[sel_bc, "scevan_label"],
        Group   = md[sel_bc, "group"],
        Cluster = as.character(md[sel_bc, "seurat_clusters"]),
        row.names = sel_bc
      )

      label_cols   <- c("tumor" = "#E41A1C", "non_tumor" = "#377EB8", "ref_normal" = "#999999")
      group_levels <- sort(unique(anno_df$Group))
      group_cols   <- scales::hue_pal()(length(group_levels))
      names(group_cols) <- group_levels
      cluster_cols <- scales::hue_pal()(length(unique(anno_df$Cluster)))
      names(cluster_cols) <- sort(unique(anno_df$Cluster))

      ha <- HeatmapAnnotation(
        Label   = anno_df$Label,
        Group   = anno_df$Group,
        Cluster = anno_df$Cluster,
        col = list(Label = label_cols, Group = group_cols, Cluster = cluster_cols),
        show_legend = TRUE,
        annotation_name_side = "left"
      )

      # Order cells: tumor first, then by cluster
      cell_order <- sel_bc[order(anno_df$Label, anno_df$Cluster)]
      cna_sub <- cna_sub[, cell_order]

      col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

      pdf(file.path(PLOT_DIR, "04_cna_heatmap.pdf"), width = 16, height = 10)
      draw(Heatmap(
        cna_sub,
        name = "CNA",
        col = col_fun,
        top_annotation = ha[cell_order, ],
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        use_raster = TRUE,
        raster_quality = 3,
        column_title = sprintf("SCEVAN CNA Profile (%d genes × %d cells)", nrow(cna_sub), ncol(cna_sub)),
        row_title = "Genomic position (ordered genes)"
      ))
      dev.off()

      png(file.path(PLOT_DIR, "04_cna_heatmap.png"), width = 1600, height = 1000, res = 150)
      draw(Heatmap(
        cna_sub,
        name = "CNA",
        col = col_fun,
        top_annotation = ha[cell_order, ],
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        use_raster = TRUE,
        raster_quality = 3,
        column_title = sprintf("SCEVAN CNA Profile (%d genes × %d cells)", nrow(cna_sub), ncol(cna_sub)),
        row_title = "Genomic position (ordered genes)"
      ))
      dev.off()

      cat("  ✅ 04_cna_heatmap\n")
    }
  } else {
    cat("  ⚠️  CNA_mtx_relat not found in .RData — skipping heatmap\n")
  }

  rm(env); gc()
} else {
  cat("  ⚠️  No .RData file found — skipping CNA heatmap\n")
}

# ============================================================================
# 7. Summary
# ============================================================================
gc()

n_tumor    <- sum(md$scevan_label == "tumor")
n_nontumor <- sum(md$scevan_label != "tumor")

cat("\n")
cat("╔══════════════════════════════════════════════════════════╗\n")
cat("║                  Validation Summary                      ║\n")
cat("╠══════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Total cells:       %6d                                ║\n", ncol(sobj)))
cat(sprintf("║  Tumor:             %6d  (%.1f%%)                       ║\n",
            n_tumor, n_tumor/ncol(sobj)*100))
cat(sprintf("║  Non-tumor:         %6d  (%.1f%%)                       ║\n",
            n_nontumor, n_nontumor/ncol(sobj)*100))
cat(sprintf("║  Tumor clusters:    %-40s ║\n", paste(tumor_clusters, collapse=", ")))
cat(sprintf("║  Ref clusters OK:   %-40s ║\n",
            ifelse(nrow(ref_tumor) > 0 && max(ref_tumor$pct_tumor) <= 5, "YES", "CHECK")))
cat("╠══════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Reports: %s\n", REPORT_DIR))
cat(sprintf("║  Plots:   %s\n", PLOT_DIR))
cat("╚══════════════════════════════════════════════════════════╝\n\n")
