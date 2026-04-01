#!/usr/bin/env Rscript
# ============================================================================
# 01_soupx_ambient_removal.R — Per-sample ambient RNA removal using SoupX
#
# Input:  Cell Ranger raw + filtered matrices + clusters (per sample)
# Output: Corrected count matrices (.rds) → results/scrna/02_qc/soupx/
#
# Usage:
#   Rscript workflow/scrna/02_qc/01_soupx_ambient_removal.R
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(SoupX)
  library(Matrix)
  library(ggplot2)
  library(data.table)
})

# Layer 2 (auto-loads Layer 1)
source(here("workflow", "scrna", "functions", "io_scrna.R"))
source(here("workflow", "scrna", "functions", "qc_utils.R"))

# ---- Config ----
SAMPLES    <- load_sample_list()
QC_PARAMS  <- load_qc_params()
out_dirs   <- scrna_output_dirs("02_qc")
soupx_plot_dir <- file.path(out_dirs$plots, "soupx")
dir.create(soupx_plot_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n╔══════════════════════════════════════════╗\n")
cat("║   Step 01: SoupX Ambient RNA Removal     ║\n")
cat(sprintf("║   Samples: %d                             ║\n", length(SAMPLES)))
cat(sprintf("║   Species: %-28s║\n", QC_PARAMS$species))
cat("╚══════════════════════════════════════════╝\n\n")

# ---- Per-sample processing ----
soupx_summary <- list()

for (sid in SAMPLES) {
  log_msg(sprintf("========== SoupX: %s ==========", sid))
  
  tryCatch({
    # Load matrices
    filt_mat <- load_filtered_matrix(sid)
    raw_mat  <- load_raw_matrix(sid)
    clusters <- load_clusters(sid)
    
    # Ensure matching genes
    common_genes <- intersect(rownames(raw_mat), rownames(filt_mat))
    raw_mat  <- raw_mat[common_genes, ]
    filt_mat <- filt_mat[common_genes, ]
    
    # Create SoupChannel
    sc <- SoupChannel(raw_mat, filt_mat)
    
    # Add clustering
    if (!is.null(clusters)) {
      common_bc <- intersect(colnames(filt_mat), names(clusters))
      if (length(common_bc) > 0.9 * ncol(filt_mat)) {
        sc <- setClusters(sc, clusters[common_bc])
        log_msg(sprintf("%s: %d/%d cells matched to clusters",
                        sid, length(common_bc), ncol(filt_mat)))
      }
    }
    
    # Estimate contamination
    sc <- autoEstCont(sc, verbose = TRUE)
    contam_frac <- sc$fit$rhoEst
    log_msg(sprintf("%s: contamination = %.3f (%.1f%%)", sid, contam_frac, contam_frac * 100))
    
    # Correct
    corrected <- adjustCounts(sc)
    
    # Save
    out_file <- file.path(out_dirs$soupx, sprintf("%s_soupx_corrected.rds", sid))
    saveRDS(corrected, out_file)
    log_msg(sprintf("%s: saved → %s", sid, basename(out_file)))
    
    # Diagnostic plots
    pdf(file.path(soupx_plot_dir, sprintf("%s_soupx_diagnostics.pdf", sid)),
        width = 12, height = 8)
    tryCatch({
      plot(sc, "rhoEst")
      title(main = sprintf("%s — Contamination Estimate (rho=%.3f)", sid, contam_frac))
    }, error = function(e) { plot.new(); text(0.5, 0.5, paste("Plot failed:", e$message)) })
    tryCatch({
      plotChangeMap(sc, corrected, geneSet = head(rownames(sc$soupProfile), 20))
      title(main = sprintf("%s — Top 20 Soup Genes", sid))
    }, error = function(e) { plot.new(); text(0.5, 0.5, paste("Plot failed:", e$message)) })
    dev.off()
    
    # Summary stats
    total_orig <- sum(filt_mat)
    total_corr <- sum(corrected)
    pct_removed <- round((1 - total_corr / total_orig) * 100, 2)
    
    gene_diff <- rowSums(filt_mat[, colnames(corrected)]) - rowSums(corrected)
    top_soup <- head(sort(gene_diff, decreasing = TRUE), 5)
    
    soupx_summary[[sid]] <- data.frame(
      sample            = sid,
      contamination_est = round(contam_frac, 4),
      total_umi_before  = total_orig,
      total_umi_after   = total_corr,
      pct_umi_removed   = pct_removed,
      n_cells           = ncol(corrected),
      n_genes           = nrow(corrected),
      top_soup_gene_1   = names(top_soup)[1],
      top_soup_gene_2   = names(top_soup)[2],
      top_soup_gene_3   = names(top_soup)[3],
      stringsAsFactors   = FALSE
    )
    
    log_msg(sprintf("%s: ✅ done (%.1f%% UMIs removed)", sid, pct_removed))
    
  }, error = function(e) {
    log_msg(sprintf("%s: ❌ FAILED — %s", sid, e$message), "ERROR")
    soupx_summary[[sid]] <<- data.frame(
      sample = sid, contamination_est = NA, total_umi_before = NA,
      total_umi_after = NA, pct_umi_removed = NA, n_cells = NA,
      n_genes = NA, top_soup_gene_1 = NA, top_soup_gene_2 = NA,
      top_soup_gene_3 = NA, stringsAsFactors = FALSE
    )
  })
  
  gc()
}

# ---- Save summary ----
summary_df <- do.call(rbind, soupx_summary)
fwrite(summary_df, file.path(out_dirs$reports, "soupx_summary.csv"))
log_msg(sprintf("Summary → %s", file.path(out_dirs$reports, "soupx_summary.csv")))

cat("\n")
print(summary_df)
cat("\n========== SoupX complete ==========\n")
