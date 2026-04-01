#!/usr/bin/env Rscript
# ============================================================================
# 02_qc_filter.R — MAD-based per-sample QC filtering
#
# Input:  Raw Seurat objects (from step 01)
# Output: QC-filtered Seurat objects → results/scrna/02_qc/qc_filtered/
#
# Method: Median Absolute Deviation (MAD)
#   nFeature_RNA:  log10 space, both directions, n_mad=3
#   nCount_RNA:    log10 space, both directions, n_mad=3
#   percent.mt:    linear space, upper only,     n_mad=3
#   Hard floors:   min_genes=200, min_umi=500, max_mt=20%
#
# Doublets NOT removed here — that's Step 03.
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(data.table)
})

source(here("workflow", "scrna", "functions", "io_scrna.R"))
source(here("workflow", "scrna", "functions", "qc_utils.R"))

# ---- Config ----
SAMPLES   <- load_sample_list()
QC_PARAMS <- load_qc_params()
out_dirs  <- scrna_output_dirs("02_qc")

raw_dir <- file.path(out_dirs$base, "raw_seurat")
qc_dir  <- file.path(out_dirs$base, "qc_filtered")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n╔══════════════════════════════════════════════════╗\n")
cat("║   Step 02: MAD-based QC Filtering                ║\n")
cat(sprintf("║   n_mad: %d   min_genes: %d   max_mt: %.0f%%       ║\n",
            QC_PARAMS$n_mad, QC_PARAMS$min_genes, QC_PARAMS$max_mt_pct))
cat("╚══════════════════════════════════════════════════╝\n\n")

filter_summary <- list()
threshold_reports <- list()

for (sid in SAMPLES) {
  log_msg(sprintf("========== MAD QC: %s ==========", sid))

  tryCatch({
    # Load raw Seurat
    sobj <- readRDS(file.path(raw_dir, sprintf("%s_raw.rds", sid)))
    n_before <- ncol(sobj)

    # Run MAD-based outlier detection
    qc_result <- mad_qc_filter(sobj, QC_PARAMS)

    # Store threshold report with sample label
    report <- qc_result$report
    report$sample <- sid
    threshold_reports[[sid]] <- report

    # Log thresholds
    for (i in seq_len(nrow(report))) {
      r <- report[i, ]
      if (r$metric == "percent.mt") {
        log_msg(sprintf("  %s: upper=%.1f (median=%.1f, n_out=%d)",
                        r$metric, r$upper_thresh, r$median_val, r$n_outlier))
      } else {
        log_msg(sprintf("  %s: [%.0f, %.0f] (median=%.0f, n_out=%d)",
                        r$metric, r$lower_thresh, r$upper_thresh,
                        r$median_val, r$n_outlier))
      }
    }

    # Apply filter
    sobj_filtered <- subset(sobj, cells = colnames(sobj)[qc_result$pass])
    n_after <- ncol(sobj_filtered)

    # Breakdown per metric (individual contribution)
    md <- sobj@meta.data
    nf_out <- detect_outliers_mad(md$nFeature_RNA, QC_PARAMS$n_mad,
                                   log_transform = TRUE, direction = "both",
                                   hard_min = QC_PARAMS$min_genes)
    nc_out <- detect_outliers_mad(md$nCount_RNA, QC_PARAMS$n_mad,
                                   log_transform = TRUE, direction = "both",
                                   hard_min = QC_PARAMS$min_umi)
    mt_out <- detect_outliers_mad(md$percent.mt, QC_PARAMS$n_mad,
                                   log_transform = FALSE, direction = "upper",
                                   hard_max = QC_PARAMS$max_mt_pct)

    log_msg(sprintf("%s: %d → %d cells (removed %d = %.1f%%)",
                    sid, n_before, n_after, n_before - n_after,
                    (n_before - n_after) / n_before * 100))

    # Save both versions
    saveRDS(sobj_filtered, file.path(qc_dir, sprintf("%s_qc_filtered.rds", sid)))
    # Keep raw version (already saved, but copy reference for viz)

    # Summary
    post_md <- sobj_filtered@meta.data
    filter_summary[[sid]] <- data.frame(
      sample            = sid,
      group             = unname(sobj$group[1]),
      n_before          = n_before,
      n_after           = n_after,
      pct_removed       = round((n_before - n_after) / n_before * 100, 2),
      removed_nFeature  = sum(nf_out),
      removed_nCount    = sum(nc_out),
      removed_mt        = sum(mt_out),
      nFeature_lower    = round(attr(nf_out, "lower")),
      nFeature_upper    = round(attr(nf_out, "upper")),
      nCount_lower      = round(attr(nc_out, "lower")),
      nCount_upper      = round(attr(nc_out, "upper")),
      mt_upper          = round(attr(mt_out, "upper"), 2),
      post_median_genes = median(post_md$nFeature_RNA),
      post_median_umi   = median(post_md$nCount_RNA),
      post_median_mt    = round(median(post_md$percent.mt), 2),
      stringsAsFactors   = FALSE
    )

    log_msg(sprintf("%s: ✅ done", sid))

  }, error = function(e) {
    log_msg(sprintf("%s: ❌ FAILED — %s", sid, e$message), "ERROR")
    filter_summary[[sid]] <<- data.frame(
      sample = sid, group = NA, n_before = NA, n_after = NA,
      pct_removed = NA, removed_nFeature = NA, removed_nCount = NA,
      removed_mt = NA, nFeature_lower = NA, nFeature_upper = NA,
      nCount_lower = NA, nCount_upper = NA, mt_upper = NA,
      post_median_genes = NA, post_median_umi = NA, post_median_mt = NA,
      stringsAsFactors = FALSE
    )
  })

  gc()
}

# Save summaries
summary_df <- do.call(rbind, filter_summary)
fwrite(summary_df, file.path(out_dirs$reports, "qc_filter_summary.csv"))
log_msg(sprintf("Filter summary → %s", file.path(out_dirs$reports, "qc_filter_summary.csv")))

# Save threshold details
thresh_df <- do.call(rbind, threshold_reports)
fwrite(thresh_df, file.path(out_dirs$reports, "mad_thresholds.csv"))
log_msg(sprintf("MAD thresholds → %s", file.path(out_dirs$reports, "mad_thresholds.csv")))

cat("\n=== Per-sample MAD Thresholds ===\n")
print(thresh_df)
cat("\n=== Filter Summary ===\n")
print(summary_df)
cat("\n========== MAD QC filtering complete ==========\n")
