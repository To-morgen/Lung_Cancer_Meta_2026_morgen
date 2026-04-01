#!/usr/bin/env Rscript
# ============================================================================
# 03_doublet_removal.R — Per-sample doublet detection on QC-filtered cells
#
# Input:  QC-filtered Seurat objects (from step 02)
# Output: Final clean Seurat objects → results/scrna/02_qc/clean/
#
# Runs scDblFinder on CLEAN cells (after QC filtering).
# This is more accurate than running on unfiltered data because:
#   - Dead/dying cells won't confuse doublet simulation
#   - Doublet rate estimation is based on true cell count
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(Seurat)
  library(Matrix)
  library(data.table)
})

source(here("workflow", "scrna", "functions", "io_scrna.R"))
source(here("workflow", "scrna", "functions", "qc_utils.R"))

# ---- Config ----
SAMPLES   <- load_sample_list()
QC_PARAMS <- load_qc_params()
out_dirs  <- scrna_output_dirs("02_qc")

qc_dir    <- file.path(out_dirs$base, "qc_filtered")
clean_dir <- file.path(out_dirs$base, "clean")
dir.create(clean_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(QC_PARAMS$doublet_seed)

cat("\n╔══════════════════════════════════════════╗\n")
cat("║   Step 03: scDblFinder Doublet Removal    ║\n")
cat(sprintf("║   Samples: %d    Seed: %d                ║\n", length(SAMPLES), QC_PARAMS$doublet_seed))
cat("╚══════════════════════════════════════════╝\n\n")

doublet_summary <- list()

for (sid in SAMPLES) {
  log_msg(sprintf("========== scDblFinder: %s ==========", sid))

  tryCatch({
    # Load QC-filtered Seurat object
    sobj <- readRDS(file.path(qc_dir, sprintf("%s_qc_filtered.rds", sid)))
    n_before <- ncol(sobj)

    log_msg(sprintf("%s: %d QC-filtered cells input", sid, n_before))

    # Extract counts for SCE
    mat <- GetAssayData(sobj, layer = "counts")

    # Create SCE
    sce <- SingleCellExperiment(assays = list(counts = mat))

    # Run scDblFinder
    sce <- scDblFinder(sce, verbose = TRUE)

    # Extract results
    labels <- sce$scDblFinder.class
    scores <- sce$scDblFinder.score

    n_singlet <- sum(labels == "singlet")
    n_doublet <- sum(labels == "doublet")
    pct_dbl   <- round(n_doublet / ncol(sce) * 100, 2)

    log_msg(sprintf("%s: %d singlets, %d doublets (%.1f%%)",
                    sid, n_singlet, n_doublet, pct_dbl))

    # Add to Seurat metadata
    sobj$scDblFinder_class <- as.character(labels)
    sobj$scDblFinder_score <- scores

    # Save doublet calls separately (for later reference)
    dbl_df <- data.frame(
      barcode           = colnames(sce),
      scDblFinder_class = as.character(labels),
      scDblFinder_score = round(scores, 4),
      stringsAsFactors  = FALSE
    )
    saveRDS(dbl_df, file.path(out_dirs$doublets, sprintf("%s_doublet_calls.rds", sid)))

    # Remove doublets → final clean object
    sobj_clean <- subset(sobj, scDblFinder_class == "singlet")
    n_after <- ncol(sobj_clean)

    log_msg(sprintf("%s: %d → %d cells after doublet removal", sid, n_before, n_after))

    # Save clean object
    saveRDS(sobj_clean, file.path(clean_dir, sprintf("%s_clean.rds", sid)))

    # Save pre-doublet-removal version (with labels) for visualization
    saveRDS(sobj, file.path(qc_dir, sprintf("%s_qc_with_doublet_labels.rds", sid)))

    doublet_summary[[sid]] <- data.frame(
      sample       = sid,
      group        = unname(sobj$group[1]),
      n_input      = n_before,
      n_singlet    = n_singlet,
      n_doublet    = n_doublet,
      pct_doublet  = pct_dbl,
      n_final      = n_after,
      median_score = round(median(scores), 4),
      stringsAsFactors = FALSE
    )

    log_msg(sprintf("%s: ✅ done", sid))

  }, error = function(e) {
    log_msg(sprintf("%s: ❌ FAILED — %s", sid, e$message), "ERROR")
    doublet_summary[[sid]] <<- data.frame(
      sample = sid, group = NA, n_input = NA, n_singlet = NA,
      n_doublet = NA, pct_doublet = NA, n_final = NA, median_score = NA,
      stringsAsFactors = FALSE
    )
  })

  gc()
}

# Save summary
summary_df <- do.call(rbind, doublet_summary)
fwrite(summary_df, file.path(out_dirs$reports, "doublet_summary.csv"))
log_msg(sprintf("Summary → %s", file.path(out_dirs$reports, "doublet_summary.csv")))

cat("\n"); print(summary_df)
cat("\n========== Doublet removal complete ==========\n")
