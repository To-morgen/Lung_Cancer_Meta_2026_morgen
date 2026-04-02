#!/usr/bin/env Rscript
# ============================================================================
# 03_sctransform.R — SCTransform v2 normalization
#
# Input:  merged_cc_scored.rds (merged + CC scores)
# Output: merged_sct.rds (SCTransform-normalized)
#
# Design decisions:
#   - vst.flavor = "v2" (glmGamPoi, more robust)
#   - vars.to.regress: default NULL (no regression)
#     → After Step 05 assessment, if CC is a problem, re-run with:
#       vars.to.regress = c("S.Score", "G2M.Score")
#       or c("CC.Difference") to preserve differentiation signal
#   - n_variable_features = 3000 (generous for heterogeneous tumor data)
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(glmGamPoi)   # required for vst.flavor = "v2"
  library(data.table)
})

source(here("workflow", "scrna", "functions", "io_scrna.R"))
source(here("workflow", "scrna", "functions", "qc_utils.R"))

# ---- Future memory limit (SCTransform on large datasets) ----
options(future.globals.maxSize = 10 * 1024^3)   # 10 GiB

# ---- Config ----
QC_PARAMS <- load_qc_params()
out_dirs  <- scrna_output_dirs("03_normalize")

# SCTransform parameters
N_FEATURES   <- QC_PARAMS$normalization$n_variable_features %||% 3000
# vars_to_regress: YAML [] → list() in R; treat empty list same as NULL
VARS_REGRESS <- QC_PARAMS$normalization$vars_to_regress
if (is.null(VARS_REGRESS) || length(VARS_REGRESS) == 0) VARS_REGRESS <- NULL
SCT_FLAVOR   <- "v2"

cat("\n╔══════════════════════════════════════════╗\n")
cat("║   Step 03: SCTransform v2                ║\n")
cat(sprintf("║   Features: %d   Regress: %-14s║\n",
    N_FEATURES, ifelse(is.null(VARS_REGRESS), "none", paste(VARS_REGRESS, collapse = ","))))
cat("╚══════════════════════════════════════════╝\n\n")

# ---- Load ----
in_file <- file.path(out_dirs$merged, "merged_cc_scored.rds")
if (!file.exists(in_file)) stop("Input not found: ", in_file)
sobj <- readRDS(in_file)
log_msg(sprintf("Loaded: %d cells × %d genes", ncol(sobj), nrow(sobj)))

# ---- SCTransform ----
log_msg("Running SCTransform v2 (this may take 10-20 min)...")
t0 <- Sys.time()

sobj <- SCTransform(
  sobj,
  vst.flavor          = SCT_FLAVOR,
  variable.features.n = N_FEATURES,
  vars.to.regress     = VARS_REGRESS,
  return.only.var.genes = FALSE,
  verbose             = TRUE,
  seed.use            = 42
)

elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
log_msg(sprintf("SCTransform complete in %s min", elapsed))
log_msg(sprintf("Default assay: %s", DefaultAssay(sobj)))
log_msg(sprintf("Variable features: %d", length(VariableFeatures(sobj)))  )

# ---- Save ----
out_file <- file.path(out_dirs$normalized, "merged_sct.rds")
saveRDS(sobj, out_file)
log_msg(sprintf("Saved → %s", out_file))

# ---- Summary ----
summary_info <- data.frame(
  metric = c("n_cells", "n_genes", "n_variable_features",
             "sct_flavor", "vars_regressed", "elapsed_min"),
  value  = c(ncol(sobj), nrow(sobj), length(VariableFeatures(sobj)),
             SCT_FLAVOR,
             ifelse(is.null(VARS_REGRESS), "none", paste(VARS_REGRESS, collapse = ",")),
             as.character(elapsed))
)
fwrite(summary_info, file.path(out_dirs$reports, "sctransform_summary.csv"))
cat("\n"); print(summary_info); cat("\n")

cat("\n========== SCTransform complete ==========\n")
gc()
