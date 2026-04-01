#!/usr/bin/env Rscript
# ============================================================================
# 01_create_seurat.R — Create per-sample Seurat objects from SoupX output
#
# Input:  SoupX-corrected matrices
# Output: Raw Seurat objects with QC metrics → results/scrna/02_qc/raw_seurat/
#
# QC metrics added: nFeature_RNA, nCount_RNA, percent.mt, percent.ribo, percent.hb
# No filtering applied here — just object creation + metric computation.
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(Matrix)
  library(data.table)
})

source(here("workflow", "scrna", "functions", "io_scrna.R"))
source(here("workflow", "scrna", "functions", "qc_utils.R"))

# ---- Config ----
SAMPLES       <- load_sample_list()
SAMPLE_GROUPS <- load_sample_groups()
QC_PARAMS     <- load_qc_params()
out_dirs      <- scrna_output_dirs("02_qc")

# Additional output dir for raw (unfiltered) Seurat objects
raw_dir <- file.path(out_dirs$base, "raw_seurat")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n╔══════════════════════════════════════════╗\n")
cat("║   Step 01: Create Seurat + QC Metrics    ║\n")
cat(sprintf("║   Samples: %d   Species: %-16s║\n", length(SAMPLES), QC_PARAMS$species))
cat("╚══════════════════════════════════════════╝\n\n")

create_summary <- list()

for (sid in SAMPLES) {
  log_msg(sprintf("========== Creating Seurat: %s ==========", sid))

  tryCatch({
    # Load SoupX-corrected matrix
    soupx_file <- file.path(out_dirs$soupx, sprintf("%s_soupx_corrected.rds", sid))
    if (!file.exists(soupx_file)) stop("SoupX output not found: ", soupx_file)
    mat <- readRDS(soupx_file)

    # Create Seurat object (minimal gene filter, NO cell filter)
    sobj <- CreateSeuratObject(
      counts       = mat,
      project      = sid,
      min.cells    = QC_PARAMS$min_cells_per_gene,
      min.features = 0   # ← 不过滤任何细胞，留给 Step 02 MAD 处理
    )

    # Add sample metadata
    sobj$sample_id <- sid
    sobj$group     <- SAMPLE_GROUPS[sid]

    # Add QC metrics
    sobj <- add_mito_pct(sobj, pattern = QC_PARAMS$mito_pattern)
    sobj <- add_ribo_pct(sobj, pattern = QC_PARAMS$ribo_pattern)
    sobj <- add_hemo_pct(sobj, species = QC_PARAMS$species)

    # Log2 metrics for reference
    sobj$log10_nFeature <- log10(sobj$nFeature_RNA + 1)
    sobj$log10_nCount   <- log10(sobj$nCount_RNA + 1)

    # Save
    out_file <- file.path(raw_dir, sprintf("%s_raw.rds", sid))
    saveRDS(sobj, out_file)

    # Summary
    md <- sobj@meta.data
    create_summary[[sid]] <- data.frame(
      sample        = sid,
      group         = SAMPLE_GROUPS[sid],
      n_cells       = ncol(sobj),
      n_genes       = nrow(sobj),
      median_genes  = median(md$nFeature_RNA),
      median_umi    = median(md$nCount_RNA),
      median_mt     = round(median(md$percent.mt), 2),
      median_ribo   = round(median(md$percent.ribo), 2),
      median_hb     = round(median(md$percent.hb, na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    )

    log_msg(sprintf("%s: ✅ %d cells, %d genes, median_mt=%.1f%%",
                    sid, ncol(sobj), nrow(sobj), median(md$percent.mt)))

  }, error = function(e) {
    log_msg(sprintf("%s: ❌ FAILED — %s", sid, e$message), "ERROR")
    create_summary[[sid]] <<- data.frame(
      sample = sid, group = SAMPLE_GROUPS[sid],
      n_cells = NA, n_genes = NA, median_genes = NA,
      median_umi = NA, median_mt = NA, median_ribo = NA, median_hb = NA,
      stringsAsFactors = FALSE
    )
  })

  gc()
}

summary_df <- do.call(rbind, create_summary)
fwrite(summary_df, file.path(out_dirs$reports, "raw_seurat_summary.csv"))
cat("\n"); print(summary_df)
cat("\n========== Seurat creation complete ==========\n")
