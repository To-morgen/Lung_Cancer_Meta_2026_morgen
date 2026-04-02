#!/usr/bin/env Rscript
# ============================================================================
# 02_cell_cycle_score.R — Score cell cycle phase
#
# Approach:
#   1. Quick NormalizeData (log-norm) just for CC scoring
#   2. CellCycleScoring using Seurat's gene list (mouse-converted)
#   3. Store S.Score, G2M.Score, Phase, CC.Difference in metadata
#   4. Do NOT regress yet — assess after PCA first (Step 05)
#
# Note: SCTransform in Step 03 will re-normalize from raw counts.
#       The log-norm here is ONLY for CC scoring.
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(data.table)
})

source(here("workflow", "scrna", "functions", "io_scrna.R"))
source(here("workflow", "scrna", "functions", "qc_utils.R"))

# ---- Config ----
QC_PARAMS <- load_qc_params()
out_dirs  <- scrna_output_dirs("03_normalize")

cat("\n╔══════════════════════════════════════════╗\n")
cat("║   Step 02: Cell Cycle Scoring            ║\n")
cat(sprintf("║   Species: %-29s║\n", QC_PARAMS$species))
cat("╚══════════════════════════════════════════╝\n\n")

# ---- Load merged object ----
merged_file <- file.path(out_dirs$merged, "merged_raw.rds")
if (!file.exists(merged_file)) stop("Merged object not found: ", merged_file)
sobj <- readRDS(merged_file)
log_msg(sprintf("Loaded merged: %d cells × %d genes", ncol(sobj), nrow(sobj)))

# ---- Convert CC genes for mouse ----
# Seurat's cc.genes.updated.2019 are human symbols
# Mouse: capitalize first letter only (e.g., MCM5 → Mcm5)
s_genes_human   <- cc.genes.updated.2019$s.genes
g2m_genes_human <- cc.genes.updated.2019$g2m.genes

if (QC_PARAMS$species == "mouse") {
  s_genes   <- stringr::str_to_title(s_genes_human)
  g2m_genes <- stringr::str_to_title(g2m_genes_human)
  log_msg("Converted CC genes to mouse format (str_to_title)")
} else {
  s_genes   <- s_genes_human
  g2m_genes <- g2m_genes_human
}

# Check gene overlap
s_found   <- sum(s_genes %in% rownames(sobj))
g2m_found <- sum(g2m_genes %in% rownames(sobj))
log_msg(sprintf("S genes found: %d/%d (%.0f%%)", s_found, length(s_genes), 100 * s_found / length(s_genes)))
log_msg(sprintf("G2M genes found: %d/%d (%.0f%%)", g2m_found, length(g2m_genes), 100 * g2m_found / length(g2m_genes)))

if (s_found < 20 || g2m_found < 20) {
  warning("Low CC gene overlap! Check species setting. Proceeding anyway...")
}

# ---- Quick log-normalize for CC scoring ----
log_msg("Running NormalizeData (log-norm, for CC scoring only)...")
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

# ---- Cell Cycle Scoring ----
log_msg("Running CellCycleScoring...")
sobj <- CellCycleScoring(
  sobj,
  s.features   = s_genes,
  g2m.features = g2m_genes,
  set.ident    = FALSE
)

# Add CC.Difference (useful for optional regression)
sobj$CC.Difference <- sobj$S.Score - sobj$G2M.Score

# ---- Phase distribution ----
phase_table <- table(sobj$Phase)
phase_pct   <- round(100 * prop.table(phase_table), 1)
log_msg("Cell cycle phase distribution:")
for (ph in names(phase_table)) {
  log_msg(sprintf("  %s: %d cells (%.1f%%)", ph, phase_table[ph], phase_pct[ph]))
}

# Per-sample phase distribution
phase_by_sample <- as.data.frame.matrix(table(sobj$sample_id, sobj$Phase))
phase_by_sample$sample <- rownames(phase_by_sample)
phase_by_sample <- phase_by_sample[, c("sample", "G1", "G2M", "S")]
phase_by_sample$total <- rowSums(phase_by_sample[, c("G1", "G2M", "S")])
phase_by_sample$pct_cycling <- round(100 * (phase_by_sample$G2M + phase_by_sample$S) / phase_by_sample$total, 1)

fwrite(phase_by_sample, file.path(out_dirs$reports, "cell_cycle_distribution.csv"))
log_msg(sprintf("CC distribution → %s", file.path(out_dirs$reports, "cell_cycle_distribution.csv")))
cat("\n"); print(phase_by_sample); cat("\n")

# ---- Save (overwrite merged with CC scores added) ----
out_file <- file.path(out_dirs$merged, "merged_cc_scored.rds")
saveRDS(sobj, out_file)
log_msg(sprintf("Saved → %s (with CC scores, pre-SCTransform)", out_file))

cat("\n========== Cell Cycle Scoring complete ==========\n")
gc()
