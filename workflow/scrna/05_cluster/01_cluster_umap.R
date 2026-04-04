#!/usr/bin/env Rscript
# ============================================================================
# 01_cluster_umap.R — FindNeighbors + FindClusters + RunUMAP
#
# Input:  Harmony-corrected Seurat object
# Output: Clustered object with UMAP at multiple resolutions
#
# Uses Harmony reduction (NOT raw PCA) for neighbor graph and UMAP.
# Tests multiple clustering resolutions for downstream selection.
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

source(here("workflow", "scrna", "functions", "io_scrna.R"))
source(here("workflow", "scrna", "functions", "qc_utils.R"))
source(here("scripts", "utils", "utils_plotting.R"))

# ---- Config ----
QC_PARAMS   <- load_qc_params()
int_dirs    <- scrna_integrate_dirs()
clust_dirs  <- scrna_cluster_dirs()

N_PCS       <- QC_PARAMS$integration$n_pcs %||% 30
RESOLUTIONS <- QC_PARAMS$clustering$resolutions %||% c(0.3, 0.5, 0.8, 1.0, 1.2)
DEFAULT_RES <- QC_PARAMS$clustering$default_resolution %||% 0.8
ALGORITHM   <- QC_PARAMS$clustering$algorithm %||% 4  # 4 = Leiden
UMAP_NEIGHBORS <- QC_PARAMS$reduction$umap_n_neighbors %||% 30
UMAP_MIN_DIST  <- QC_PARAMS$reduction$umap_min_dist %||% 0.3

cat("\n╔══════════════════════════════════════════════════════╗\n")
cat("║   Phase 5 Step 01: Clustering + UMAP                 ║\n")
cat(sprintf("║   Dims: 1:%d   Algorithm: %s                       ║\n",
            N_PCS, ifelse(ALGORITHM == 4, "Leiden", "Louvain")))
cat(sprintf("║   Resolutions: %s             ║\n",
            paste(RESOLUTIONS, collapse = ", ")))
cat("╚══════════════════════════════════════════════════════╝\n\n")

# ---- Load Harmony object ----
harmony_file <- file.path(int_dirs$harmony, "seurat_harmony.rds")
if (!file.exists(harmony_file)) stop("Harmony object not found: ", harmony_file)

log_msg(sprintf("Loading: %s", basename(harmony_file)))
sobj <- readRDS(harmony_file)
log_msg(sprintf("  Cells: %d, Reductions: %s",
                ncol(sobj), paste(names(sobj@reductions), collapse = ", ")))

stopifnot("harmony" %in% names(sobj@reductions))

# ---- FindNeighbors ----
log_msg(sprintf("FindNeighbors on Harmony 1:%d ...", N_PCS))
t0 <- Sys.time()

sobj <- FindNeighbors(
  sobj,
  reduction = "harmony",
  dims      = 1:N_PCS,
  verbose   = TRUE
)

log_msg(sprintf("  FindNeighbors done in %.1f min",
                difftime(Sys.time(), t0, units = "mins")))

# ---- FindClusters at multiple resolutions ----
log_msg(sprintf("FindClusters at %d resolutions...", length(RESOLUTIONS)))

for (res in RESOLUTIONS) {
  log_msg(sprintf("  Resolution %.1f ...", res))
  sobj <- FindClusters(
    sobj,
    resolution = res,
    algorithm  = ALGORITHM,
    verbose    = FALSE
  )
  # Seurat auto-names: SCT_snn_res.0.3, etc.
  # Also create a clean column name
  col_name <- sprintf("clusters_res%.1f", res)
  auto_col <- grep(sprintf("res\\.%s$", res), colnames(sobj@meta.data), value = TRUE)
  if (length(auto_col) > 0) {
    sobj@meta.data[[col_name]] <- sobj@meta.data[[auto_col[1]]]
  }
  n_clusters <- length(unique(sobj@meta.data[[col_name]]))
  log_msg(sprintf("    res=%.1f → %d clusters", res, n_clusters))
}

# Set default active ident to default resolution
default_col <- sprintf("clusters_res%.1f", DEFAULT_RES)
if (default_col %in% colnames(sobj@meta.data)) {
  Idents(sobj) <- default_col
  log_msg(sprintf("Active Idents set to %s (%d clusters)",
                  default_col, length(levels(Idents(sobj)))))
}

# ---- RunUMAP ----
log_msg(sprintf("RunUMAP on Harmony 1:%d (n_neighbors=%d, min_dist=%.2f)...",
                N_PCS, UMAP_NEIGHBORS, UMAP_MIN_DIST))
t0 <- Sys.time()

sobj <- RunUMAP(
  sobj,
  reduction   = "harmony",
  dims        = 1:N_PCS,
  n.neighbors = UMAP_NEIGHBORS,
  min.dist    = UMAP_MIN_DIST,
  verbose     = TRUE
)

log_msg(sprintf("  RunUMAP done in %.1f min",
                difftime(Sys.time(), t0, units = "mins")))

# ---- Basic UMAP plots ----
log_msg("Generating UMAP plots...")

# By cluster (default resolution)
p_clust <- DimPlot(sobj, reduction = "umap", label = TRUE, label.size = 3, pt.size = 0.1) +
  labs(title = sprintf("UMAP — Clusters (res=%.1f, %s)",
                       DEFAULT_RES, ifelse(ALGORITHM == 4, "Leiden", "Louvain"))) +
  theme_project() + NoLegend()

# By sample
p_sample <- DimPlot(sobj, reduction = "umap", group.by = "sample_id", pt.size = 0.1) +
  labs(title = "UMAP — by Sample") +
  theme_project()

# By group
p_group <- DimPlot(sobj, reduction = "umap", group.by = "group", pt.size = 0.1) +
  labs(title = "UMAP — by Group") +
  theme_project()

# By cell cycle phase
if ("Phase" %in% colnames(sobj@meta.data)) {
  p_phase <- DimPlot(sobj, reduction = "umap", group.by = "Phase", pt.size = 0.1) +
    labs(title = "UMAP — by Cell Cycle Phase") +
    theme_project()
} else {
  p_phase <- ggplot() + theme_void() + labs(title = "Phase not available")
}

# Save combined
pdf(file.path(clust_dirs$plots, "01_umap_overview.pdf"), width = 16, height = 14)
print((p_clust | p_sample) / (p_group | p_phase))
dev.off()
png(file.path(clust_dirs$plots, "01_umap_overview.png"), width = 1600, height = 1400, res = 150)
print((p_clust | p_sample) / (p_group | p_phase))
dev.off()
log_msg("  ✅ 01_umap_overview saved")

# ---- Multi-resolution comparison ----
log_msg("Plotting multi-resolution comparison...")

res_plots <- list()
for (res in RESOLUTIONS) {
  col <- sprintf("clusters_res%.1f", res)
  if (col %in% colnames(sobj@meta.data)) {
    n_cl <- length(unique(sobj@meta.data[[col]]))
    res_plots[[sprintf("res_%.1f", res)]] <-
      DimPlot(sobj, reduction = "umap", group.by = col,
              label = TRUE, label.size = 2.5, pt.size = 0.05) +
      labs(title = sprintf("res=%.1f (%d clusters)", res, n_cl)) +
      theme_project() + NoLegend()
  }
}

pdf(file.path(clust_dirs$plots, "02_multi_resolution.pdf"), width = 20, height = 8)
print(wrap_plots(res_plots, nrow = 1) +
        plot_annotation(title = "Clustering Resolution Comparison",
                        theme = theme(plot.title = element_text(face = "bold", size = 14))))
dev.off()
png(file.path(clust_dirs$plots, "02_multi_resolution.png"), width = 2000, height = 800, res = 150)
print(wrap_plots(res_plots, nrow = 1) +
        plot_annotation(title = "Clustering Resolution Comparison",
                        theme = theme(plot.title = element_text(face = "bold", size = 14))))
dev.off()
log_msg("  ✅ 02_multi_resolution saved")

# ---- Split by sample ----
log_msg("Plotting UMAP split by sample...")

pdf(file.path(clust_dirs$plots, "03_umap_split_sample.pdf"), width = 20, height = 10)
print(DimPlot(sobj, reduction = "umap", split.by = "sample_id",
              label = TRUE, label.size = 2, pt.size = 0.05, ncol = 3) +
        labs(title = sprintf("UMAP split by sample (res=%.1f)", DEFAULT_RES)) +
        theme_project() + NoLegend())
dev.off()
png(file.path(clust_dirs$plots, "03_umap_split_sample.png"), width = 2000, height = 1000, res = 150)
print(DimPlot(sobj, reduction = "umap", split.by = "sample_id",
              label = TRUE, label.size = 2, pt.size = 0.05, ncol = 3) +
        labs(title = sprintf("UMAP split by sample (res=%.1f)", DEFAULT_RES)) +
        theme_project() + NoLegend())
dev.off()
log_msg("  ✅ 03_umap_split_sample saved")

# ---- Save ----
out_file <- file.path(clust_dirs$objects, "seurat_clustered.rds")
log_msg(sprintf("Saving clustered object → %s", out_file))
saveRDS(sobj, out_file)
log_msg(sprintf("  File size: %.1f GB", file.info(out_file)$size / 1e9))

# ---- Summary table ----
res_summary <- data.frame(
  resolution = RESOLUTIONS,
  n_clusters = sapply(RESOLUTIONS, function(r) {
    col <- sprintf("clusters_res%.1f", r)
    if (col %in% colnames(sobj@meta.data)) length(unique(sobj@meta.data[[col]])) else NA
  }),
  stringsAsFactors = FALSE
)

fwrite(res_summary, file.path(clust_dirs$reports, "resolution_summary.csv"))

cat("\n=== Resolution vs Number of Clusters ===\n")
print(res_summary)

cat("\n========== Clustering + UMAP complete ==========\n")
