#!/usr/bin/env Rscript
# ============================================================================
# 01_harmony.R — Batch correction using Harmony
#
# Input:  Merged SCTransform + PCA object (from Phase 3)
# Output: Seurat object with "harmony" reduction
#
# Harmony operates in PCA space → does NOT alter gene expression matrix.
# Downstream (FindNeighbors, RunUMAP) should use reduction = "harmony".
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

source(here("workflow", "scrna", "functions", "io_scrna.R"))
source(here("workflow", "scrna", "functions", "qc_utils.R"))
source(here("scripts", "utils", "utils_plotting.R"))

# ---- Config ----
QC_PARAMS <- load_qc_params()
int_dirs  <- scrna_integrate_dirs()

# Key parameters
GROUP_VAR <- QC_PARAMS$integration$group_by %||% "sample_id"
N_PCS     <- QC_PARAMS$integration$n_pcs %||% 30

cat("\n╔══════════════════════════════════════════════════╗\n")
cat("║   Phase 4 Step 01: Harmony Batch Correction      ║\n")
cat(sprintf("║   Batch variable: %-30s  ║\n", GROUP_VAR))
cat(sprintf("║   PCs: 1:%d                                      ║\n", N_PCS))
cat("╚══════════════════════════════════════════════════╝\n\n")

# ---- Load PCA object ----
pca_file <- scrna_base("03_normalize", "pca", "merged_pca.rds")
if (!file.exists(pca_file)) {
  stop("PCA object not found: ", pca_file, "\nRun Phase 3 (04_pca.R) first!")
}

log_msg(sprintf("Loading PCA object: %s", basename(pca_file)))
log_msg(sprintf("  File size: %.1f GB", file.info(pca_file)$size / 1e9))

sobj <- readRDS(pca_file)
log_msg(sprintf("  Cells: %d, Genes: %d", ncol(sobj), nrow(sobj)))
log_msg(sprintf("  Batch levels (%s): %s", GROUP_VAR,
                paste(unique(sobj@meta.data[[GROUP_VAR]]), collapse = ", ")))

# Verify PCA exists
stopifnot("pca" %in% names(sobj@reductions))
log_msg(sprintf("  PCA dims available: %d", ncol(Embeddings(sobj, "pca"))))

# ---- Pre-Harmony visualization ----
log_msg("Plotting pre-Harmony PCA...")

p_pre1 <- DimPlot(sobj, reduction = "pca", group.by = GROUP_VAR, pt.size = 0.1) +
  labs(title = sprintf("PCA (pre-Harmony) — colored by %s", GROUP_VAR)) +
  theme_project()

p_pre2 <- DimPlot(sobj, reduction = "pca", group.by = "group", pt.size = 0.1) +
  labs(title = "PCA (pre-Harmony) — colored by group") +
  theme_project()

pdf(file.path(int_dirs$plots, "01_pca_pre_harmony.pdf"), width = 16, height = 7)
print(p_pre1 | p_pre2)
dev.off()
png(file.path(int_dirs$plots, "01_pca_pre_harmony.png"), width = 1600, height = 700, res = 150)
print(p_pre1 | p_pre2)
dev.off()
log_msg("  ✅ Pre-Harmony PCA plots saved")

# ---- Run Harmony ----
log_msg("Running Harmony...")
t0 <- Sys.time()

sobj <- RunHarmony(
  sobj,
  group.by.vars = GROUP_VAR,
  reduction      = "pca",
  dims.use       = 1:N_PCS,
  assay.use      = "SCT",
  reduction.save = "harmony",
  verbose        = TRUE
)

elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
log_msg(sprintf("Harmony completed in %.1f min", elapsed))
log_msg(sprintf("  Harmony dims: %d", ncol(Embeddings(sobj, "harmony"))))

# ---- Post-Harmony visualization ----
log_msg("Plotting post-Harmony embeddings...")

p_post1 <- DimPlot(sobj, reduction = "harmony", group.by = GROUP_VAR, pt.size = 0.1) +
  labs(title = sprintf("Harmony — colored by %s", GROUP_VAR)) +
  theme_project()

p_post2 <- DimPlot(sobj, reduction = "harmony", group.by = "group", pt.size = 0.1) +
  labs(title = "Harmony — colored by group") +
  theme_project()

pdf(file.path(int_dirs$plots, "02_harmony_embeddings.pdf"), width = 16, height = 7)
print(p_post1 | p_post2)
dev.off()
png(file.path(int_dirs$plots, "02_harmony_embeddings.png"), width = 1600, height = 700, res = 150)
print(p_post1 | p_post2)
dev.off()

# ---- Before vs After comparison ----
p_compare <- (p_pre1 + labs(title = "Before Harmony")) |
             (p_post1 + labs(title = "After Harmony"))

pdf(file.path(int_dirs$plots, "03_harmony_before_after.pdf"), width = 16, height = 7)
print(p_compare)
dev.off()
png(file.path(int_dirs$plots, "03_harmony_before_after.png"), width = 1600, height = 700, res = 150)
print(p_compare)
dev.off()
log_msg("  ✅ Harmony comparison plots saved")

# ---- Batch mixing metric (simple) ----
log_msg("Computing batch mixing summary...")

harmony_emb <- Embeddings(sobj, "harmony")[, 1:min(5, N_PCS)]
batch_labels <- sobj@meta.data[[GROUP_VAR]]

# Per-batch centroid distance in Harmony space
batch_centroids <- aggregate(harmony_emb, by = list(batch = batch_labels), FUN = mean)
rownames(batch_centroids) <- batch_centroids$batch
batch_centroids$batch <- NULL

centroid_dists <- as.matrix(dist(batch_centroids))

log_msg("  Batch centroid distances (Harmony 1:5):")
print(round(centroid_dists, 3))

# Save centroid distances
fwrite(as.data.frame(centroid_dists, check.names = FALSE),
       file.path(int_dirs$reports, "harmony_batch_centroid_distances.csv"),
       row.names = TRUE)

# ---- Save ----
out_file <- file.path(int_dirs$harmony, "seurat_harmony.rds")
log_msg(sprintf("Saving Harmony object → %s", out_file))
saveRDS(sobj, out_file)
log_msg(sprintf("  File size: %.1f GB", file.info(out_file)$size / 1e9))

# Summary
summary_df <- data.frame(
  parameter    = c("batch_variable", "n_pcs", "n_cells", "n_batches",
                    "harmony_dims", "elapsed_min"),
  value        = c(GROUP_VAR, N_PCS, ncol(sobj),
                    length(unique(batch_labels)),
                    ncol(Embeddings(sobj, "harmony")),
                    as.character(elapsed)),
  stringsAsFactors = FALSE
)
fwrite(summary_df, file.path(int_dirs$reports, "harmony_summary.csv"))

cat("\n"); print(summary_df)
cat("\n========== Harmony integration complete ==========\n")
