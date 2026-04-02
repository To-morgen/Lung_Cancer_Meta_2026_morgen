#!/usr/bin/env Rscript
# ============================================================================
# 01_merge_samples.R — Merge all clean per-sample Seurat objects
#
# Input:  results/scrna/02_qc/clean/*_clean.rds  (6 objects)
# Output: results/scrna/03_normalize/merged/merged_raw.rds
#
# In Seurat v5, merge creates separate layers per sample.
# We JoinLayers to create a single count matrix for SCTransform.
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(data.table)
})

source(here("workflow", "scrna", "functions", "io_scrna.R"))

# ---- Config ----
SAMPLES  <- load_sample_list()
qc_dirs  <- scrna_output_dirs("02_qc")
out_dirs <- scrna_output_dirs("03_normalize")

cat("\n╔══════════════════════════════════════════╗\n")
cat("║   Step 01: Merge Clean Seurat Objects    ║\n")
cat(sprintf("║   Samples: %-28d║\n", length(SAMPLES)))
cat("╚══════════════════════════════════════════╝\n\n")

# ---- Load all clean objects ----
obj_list <- list()
for (sid in SAMPLES) {
  f <- file.path(qc_dirs$clean, sprintf("%s_clean.rds", sid))
  if (!file.exists(f)) stop("Clean object not found: ", f)
  obj_list[[sid]] <- readRDS(f)
  log_msg(sprintf("Loaded %s: %d cells × %d genes", sid, ncol(obj_list[[sid]]), nrow(obj_list[[sid]])))
}

# ---- Merge ----
log_msg("Merging all samples...")
merged <- merge(
  x           = obj_list[[1]],
  y           = obj_list[-1],
  add.cell.ids = SAMPLES,     # prefix barcodes to avoid collision
  project     = "Lung_Cancer_Meta"
)

log_msg(sprintf("Merged object: %d cells × %d genes", ncol(merged), nrow(merged)))
log_msg(sprintf("Layers: %s", paste(Layers(merged), collapse = ", ")))

# ---- Join Layers (Seurat v5) ----
log_msg("Joining layers...")
merged <- JoinLayers(merged)
log_msg(sprintf("After JoinLayers: %s", paste(Layers(merged), collapse = ", ")))

# ---- Verify metadata ----
sample_table <- table(merged$sample_id)
log_msg("Cells per sample:")
print(sample_table)

group_table <- table(merged$group)
log_msg("Cells per group:")
print(group_table)

# ---- Save ----
out_file <- file.path(out_dirs$merged, "merged_raw.rds")
saveRDS(merged, out_file)
log_msg(sprintf("Saved → %s", out_file))

# ---- Summary report ----
summary_df <- data.frame(
  sample   = names(sample_table),
  n_cells  = as.integer(sample_table),
  stringsAsFactors = FALSE
)
summary_df$group <- unname(sapply(summary_df$sample, function(s) merged$group[merged$sample_id == s][1]))
summary_df <- rbind(summary_df, data.frame(
  sample = "TOTAL", n_cells = ncol(merged), group = "-"
))
fwrite(summary_df, file.path(out_dirs$reports, "merge_summary.csv"))
log_msg(sprintf("Summary → %s", file.path(out_dirs$reports, "merge_summary.csv")))

cat("\n"); print(summary_df)
cat("\n========== Merge complete ==========\n")

rm(obj_list); gc()
