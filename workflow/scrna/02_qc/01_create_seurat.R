#!/usr/bin/env Rscript
# ============================================================================
# 01_create_seurat.R — Create per-sample Seurat objects
#
# Input modes (auto-detected from dataset config):
#   1. SoupX-corrected matrices  → results/.../soupx/{sid}_soupx_corrected.rds
#   2. Raw 10x MTX directories   → config$mtx_dir/{sid}/
#   3. CellRanger filtered output → config$cellranger_out/{sid}/outs/filtered_feature_bc_matrix/
#
# Priority: SoupX output > MTX dir > CellRanger
#
# Output: Raw Seurat objects with QC metrics → results/scrna/.../02_qc/raw_seurat/
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
DS_CONFIG     <- load_dataset_config()
out_dirs      <- scrna_output_dirs("02_qc")

raw_dir <- file.path(out_dirs$base, "raw_seurat")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

# ── Determine input mode ──
# Check if SoupX was skipped (no soupx dir, or empty, or config says so)
soupx_available <- dir.exists(out_dirs$soupx) &&
  length(list.files(out_dirs$soupx, pattern = "_soupx_corrected\\.rds$")) > 0

mtx_dir <- DS_CONFIG$mtx_dir          # NULL if not defined
cr_out  <- DS_CONFIG$cellranger_out   # NULL if not defined

if (soupx_available) {
  INPUT_MODE <- "soupx"
} else if (!is.null(mtx_dir) && dir.exists(mtx_dir)) {
  INPUT_MODE <- "mtx"
} else if (!is.null(cr_out) && dir.exists(cr_out)) {
  INPUT_MODE <- "cellranger"
} else {
  stop("No valid input found. Need SoupX output, mtx_dir, or cellranger_out in dataset config.")
}

cat("\n╔══════════════════════════════════════════════════╗\n")
cat("║   Step 01: Create Seurat + QC Metrics            ║\n")
cat(sprintf("║   Samples: %d   Species: %-16s      ║\n", length(SAMPLES), QC_PARAMS$species))
cat(sprintf("║   Input:   %-38s║\n", INPUT_MODE))
cat("╚══════════════════════════════════════════════════╝\n\n")

log_msg(sprintf("Input mode: %s", INPUT_MODE))

create_summary <- list()

for (sid in SAMPLES) {
  log_msg(sprintf("========== Creating Seurat: %s ==========", sid))

  tryCatch({

    # ══════════════════════════════════════════════════
    # Load count matrix based on input mode
    # ══════════════════════════════════════════════════
    if (INPUT_MODE == "soupx") {
      # ── Mode 1: SoupX-corrected RDS ──
      soupx_file <- file.path(out_dirs$soupx, sprintf("%s_soupx_corrected.rds", sid))
      if (!file.exists(soupx_file)) stop("SoupX output not found: ", soupx_file)
      mat <- readRDS(soupx_file)
      log_msg(sprintf("%s: loaded SoupX-corrected matrix", sid))

    } else if (INPUT_MODE == "mtx") {
      # ── Mode 2: Pre-processed 10x MTX (e.g. GEO download) ──
      sample_mtx_dir <- file.path(mtx_dir, sid)
      if (!dir.exists(sample_mtx_dir)) stop("MTX dir not found: ", sample_mtx_dir)
      mat <- Read10X(data.dir = sample_mtx_dir)
      # Read10X may return a list for multi-modal data (e.g. Gene Expression + Antibody Capture)
      if (is.list(mat)) {
        log_msg(sprintf("%s: Read10X returned list with: %s — using 'Gene Expression'",
                        sid, paste(names(mat), collapse = ", ")))
        mat <- mat[["Gene Expression"]]
      }
      log_msg(sprintf("%s: loaded 10x MTX from %s", sid, sample_mtx_dir))

    } else if (INPUT_MODE == "cellranger") {
      # ── Mode 3: CellRanger filtered output ──
      cr_dir <- file.path(cr_out, sid, "outs", "filtered_feature_bc_matrix")
      if (!dir.exists(cr_dir)) stop("CellRanger dir not found: ", cr_dir)
      mat <- Read10X(data.dir = cr_dir)
      if (is.list(mat)) mat <- mat[["Gene Expression"]]
      log_msg(sprintf("%s: loaded CellRanger filtered matrix from %s", sid, cr_dir))
    }

    # ── FIX: Convert dgTMatrix → dgCMatrix (Seurat v5 compatibility) ──
    if (inherits(mat, "dgTMatrix") || inherits(mat, "TsparseMatrix")) {
      log_msg(sprintf("%s: Converting %s → dgCMatrix", sid, class(mat)[1]))
      mat <- as(mat, "CsparseMatrix")
    }

    log_msg(sprintf("%s: matrix %d genes × %d cells", sid, nrow(mat), ncol(mat)))
    log_msg(sprintf("%s: first 3 barcodes: %s", sid, paste(head(colnames(mat), 3), collapse = ", ")))

    # Create Seurat object (minimal gene filter, NO cell filter)
    sobj <- CreateSeuratObject(
      counts       = mat,
      project      = sid,
      min.cells    = QC_PARAMS$min_cells_per_gene,
      min.features = 0
    )

    log_msg(sprintf("%s: Seurat created — %d cells, %d genes", sid, ncol(sobj), nrow(sobj)))
    log_msg(sprintf("%s: first 3 cell names: %s", sid, paste(head(colnames(sobj), 3), collapse = ", ")))

    # ── Metadata ──
    sobj$sample_id <- sid
    sobj$group     <- unname(SAMPLE_GROUPS[sid])

    # Add QC metrics
    sobj <- add_mito_pct(sobj, pattern = QC_PARAMS$mito_pattern)
    sobj <- add_ribo_pct(sobj, pattern = QC_PARAMS$ribo_pattern)
    sobj <- add_hemo_pct(sobj, species = QC_PARAMS$species)

    # Log10 metrics for reference
    sobj$log10_nFeature <- log10(sobj$nFeature_RNA + 1)
    sobj$log10_nCount   <- log10(sobj$nCount_RNA + 1)

    # Save
    out_file <- file.path(raw_dir, sprintf("%s_raw.rds", sid))
    saveRDS(sobj, out_file)

    # Summary
    md <- sobj@meta.data
    create_summary[[sid]] <- data.frame(
      sample        = sid,
      group         = unname(SAMPLE_GROUPS[sid]),
      input_mode    = INPUT_MODE,
      n_cells       = ncol(sobj),
      n_genes       = nrow(sobj),
      median_genes  = median(md$nFeature_RNA),
      median_umi    = median(md$nCount_RNA),
      median_mt     = round(median(md$percent.mt), 2),
      median_ribo   = round(median(md$percent.ribo), 2),
      median_hb     = round(median(md$percent.hb, na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    )

    log_msg(sprintf("%s: ✅ %d cells, %d genes, median_mt=%.1f%%, median_ribo=%.1f%%",
                    sid, ncol(sobj), nrow(sobj),
                    median(md$percent.mt), median(md$percent.ribo)))

  }, error = function(e) {
    log_msg(sprintf("%s: ❌ FAILED — %s", sid, e$message), "ERROR")
    log_msg(sprintf("%s: Traceback: %s", sid,
                    paste(capture.output(traceback()), collapse = "\n")), "DEBUG")
    create_summary[[sid]] <<- data.frame(
      sample = sid, group = unname(SAMPLE_GROUPS[sid]),
      input_mode = INPUT_MODE,
      n_cells = NA, n_genes = NA, median_genes = NA,
      median_umi = NA, median_mt = NA, median_ribo = NA, median_hb = NA,
      stringsAsFactors = FALSE
    )
  })

  gc()
}

summary_df <- do.call(rbind, create_summary)
fwrite(summary_df, file.path(out_dirs$reports, "raw_seurat_summary.csv"))
log_msg(sprintf("Summary → %s", file.path(out_dirs$reports, "raw_seurat_summary.csv")))
cat("\n"); print(summary_df)
cat("\n========== Seurat creation complete ==========\n")
