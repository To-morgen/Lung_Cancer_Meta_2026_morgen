#!/usr/bin/env Rscript
# ============================================================================
# 04_pca.R — Run PCA + diagnostics
#
# Input:  merged_sct.rds
# Output: merged_pca.rds + ElbowPlot + variance explained table
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ggplot2)
  library(data.table)
})

source(here("workflow", "scrna", "functions", "io_scrna.R"))
source(here("workflow", "scrna", "functions", "qc_utils.R"))

# ---- Config ----
QC_PARAMS <- load_qc_params()
out_dirs  <- scrna_output_dirs("03_normalize")

N_PCS <- QC_PARAMS$normalization$n_pcs %||% 50

cat("\n╔══════════════════════════════════════════╗\n")
cat("║   Step 04: PCA                           ║\n")
cat(sprintf("║   Computing %d PCs                       ║\n", N_PCS))
cat("╚══════════════════════════════════════════╝\n\n")

# ---- Load ----
in_file <- file.path(out_dirs$normalized, "merged_sct.rds")
if (!file.exists(in_file)) stop("Input not found: ", in_file)
sobj <- readRDS(in_file)
log_msg(sprintf("Loaded: %d cells × %d genes", ncol(sobj), nrow(sobj)))

# ---- PCA ----
log_msg(sprintf("Running PCA (npcs=%d)...", N_PCS))
t0 <- Sys.time()
sobj <- RunPCA(sobj, npcs = N_PCS, verbose = FALSE, seed.use = 42)
elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
log_msg(sprintf("PCA complete in %s min", elapsed))

# ---- Variance explained ----
stdev    <- Stdev(sobj, reduction = "pca")
var_exp  <- (stdev^2) / sum(stdev^2) * 100
cum_var  <- cumsum(var_exp)

var_df <- data.frame(
  PC             = 1:length(stdev),
  stdev          = round(stdev, 4),
  pct_variance   = round(var_exp, 2),
  cumulative_pct = round(cum_var, 2)
)

fwrite(var_df, file.path(out_dirs$reports, "pca_variance.csv"))
log_msg(sprintf("Variance explained → %s", file.path(out_dirs$reports, "pca_variance.csv")))

# Key PCs summary
for (n in c(10, 20, 30, 40, 50)) {
  if (n <= nrow(var_df)) {
    log_msg(sprintf("  PC 1-%d: %.1f%% cumulative variance", n, var_df$cumulative_pct[n]))
  }
}

# ---- Elbow Plot ----
p_elbow <- ElbowPlot(sobj, ndims = N_PCS) +
  geom_vline(xintercept = 30, linetype = "dashed", color = "red", alpha = 0.5) +
  labs(title = "PCA Elbow Plot",
       subtitle = sprintf("%d cells, %d variable features", ncol(sobj), length(VariableFeatures(sobj)))) +
  theme_minimal()

ggsave(file.path(out_dirs$plots, "elbow_plot.pdf"), p_elbow, width = 8, height = 5)
ggsave(file.path(out_dirs$plots, "elbow_plot.png"), p_elbow, width = 8, height = 5, dpi = 150)
log_msg("Elbow plot saved")

# ---- Top genes per PC ----
sink(file.path(out_dirs$reports, "pca_top_genes.txt"))
print(sobj[["pca"]], dims = 1:10, nfeatures = 10)
sink()
log_msg("Top PCA genes → pca_top_genes.txt")

# ---- Save ----
out_file <- file.path(out_dirs$pca, "merged_pca.rds")
saveRDS(sobj, out_file)
log_msg(sprintf("Saved → %s", out_file))

cat("\nTop 15 PCs variance:\n")
print(var_df[1:15, ], row.names = FALSE)

cat("\n========== PCA complete ==========\n")
gc()
