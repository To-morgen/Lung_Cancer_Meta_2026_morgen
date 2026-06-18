#!/usr/bin/env Rscript
# ============================================================================
# 01_scevan.R — SCEVAN CNV inference: tumor/normal classification
#
# ONLY runs SCEVAN. Does NOT do validation or plotting.
# Output: a standardized Seurat object with scevan_label in metadata.
#
# Next step: 02_validate_scevan.R (reads the same output regardless of path)
#
# Two completion paths:
#   Path A (happy): SCEVAN completes → parse results → standardize → save
#   Path B (crash): SCEVAN crashes  → 03_recover → standardize → save
# Both paths produce IDENTICAL output format.
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(SCEVAN)
  library(data.table)
  library(yaml)
})

# ============================================================================
# 1. Paths & Config
# ============================================================================

proj_root <- Sys.getenv("LUNGMETA_ROOT")
if (proj_root == "") {
  proj_root <- normalizePath(file.path(here::here(), "..", ".."))
  cat(sprintf("LUNGMETA_ROOT not set, inferred: %s\n", proj_root))
}

module_root <- here::here()

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

cat("\n")
cat("╔══════════════════════════════════════════════════╗\n")
cat("║   CNV Module: SCEVAN Tumor/Normal Inference      ║\n")
cat(sprintf("║   Project:  %-36s║\n", basename(proj_root)))
cat(sprintf("║   Time:     %s            ║\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("╚══════════════════════════════════════════════════╝\n\n")

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else if (nzchar(Sys.getenv("CNV_CONFIG", unset = ""))) {
  Sys.getenv("CNV_CONFIG")
} else {
  file.path(module_root, "configs", "cnv_params.yaml")
}
if (!file.exists(cfg_path)) stop("Config not found: ", cfg_path)
cfg <- yaml::read_yaml(cfg_path)

path_from_root <- function(path) {
  if (grepl("^/", path)) path else file.path(proj_root, path)
}

scevan_cfg <- cfg$scevan
if (is.null(scevan_cfg$organism)) scevan_cfg$organism <- cfg$species
if (is.null(scevan_cfg$par_cores)) scevan_cfg$par_cores <- 8
scevan_cores_env <- Sys.getenv("SCEVAN_CORES", unset = "")
if (nzchar(scevan_cores_env)) {
  scevan_cfg$par_cores <- as.integer(scevan_cores_env)
  if (is.na(scevan_cfg$par_cores) || scevan_cfg$par_cores < 1) {
    stop("Invalid SCEVAN_CORES: ", scevan_cores_env)
  }
}
if (is.null(scevan_cfg$subclones)) scevan_cfg$subclones <- TRUE
if (is.null(scevan_cfg$sample_name)) scevan_cfg$sample_name <- cfg$dataset_id %||% "LLC_tumor"

cat(sprintf("Config:     %s\n", cfg_path))
cat(sprintf("Dataset:    %s\n", cfg$dataset_id %||% "unspecified"))
cat(sprintf("Species:    %s\n", cfg$species))
cat(sprintf("Organism:   %s\n", scevan_cfg$organism))
cat(sprintf("Cores:      %d\n", scevan_cfg$par_cores))
cat(sprintf("Subclones:  %s\n", ifelse(isTRUE(scevan_cfg$subclones), "yes", "no")))

local_out <- if (!is.null(cfg$output$module_scevan)) {
  path_from_root(cfg$output$module_scevan)
} else if (!is.null(cfg$output$module_results)) {
  file.path(module_root, cfg$output$module_results, "scevan")
} else if (!is.null(cfg$output$scevan_local)) {
  path_from_root(cfg$output$scevan_local)
} else {
  file.path(module_root, "results", cfg$dataset_id %||% "default", "scevan")
}

main_out <- if (!is.null(cfg$output$scevan)) {
  path_from_root(cfg$output$scevan)
} else {
  file.path(path_from_root(cfg$output$main_results), "scevan")
}

report_out <- if (!is.null(cfg$output$reports)) {
  path_from_root(cfg$output$reports)
} else {
  file.path(path_from_root(cfg$output$main_results), "reports")
}

for (d in c(local_out, main_out, report_out)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

Sys.setenv(CNV_CONFIG = normalizePath(cfg_path), CNV_LOCAL_SCEVAN_DIR = local_out)

FINAL_OUTPUT <- file.path(main_out, "seurat_with_scevan.rds")

# ============================================================================
# 2. Load Seurat Object
# ============================================================================

input_path <- path_from_root(cfg$input$seurat_object)
if (!file.exists(input_path)) stop("Input not found: ", input_path)

cat(sprintf("\nLoading: %s\n", basename(input_path)))
sobj <- readRDS(input_path)

tryCatch(
  sobj <- JoinLayers(sobj),
  error = function(e) cat(sprintf("  JoinLayers skipped: %s\n", e$message))
)

cat(sprintf("  Cells: %d | Genes: %d | Clusters: %s\n",
            ncol(sobj), nrow(sobj),
            paste(sort(unique(sobj$seurat_clusters)), collapse = ",")))

# ============================================================================
# 3. Optional Downsample
# ============================================================================

ds_n <- scevan_cfg$downsample
if (!is.null(ds_n) && is.numeric(ds_n) && ds_n < ncol(sobj)) {
  cat(sprintf("\n⚠️  Downsampling: %d → %d cells\n", ncol(sobj), ds_n))
  set.seed(42)
  sobj <- subset(sobj, cells = sample(colnames(sobj), ds_n))
}

# ============================================================================
# 4. Extract Count Matrix (SPARSE)
# ============================================================================

cat("\nExtracting count matrix (sparse)...\n")
counts <- GetAssayData(sobj, assay = "RNA", layer = "counts")
cat(sprintf("  %d genes × %d cells (%s)\n",
            nrow(counts), ncol(counts), class(counts)[1]))

# ============================================================================
# 5. Reference Cells
# ============================================================================

ref_config <- scevan_cfg$reference_clusters

if (is.null(ref_config) || (is.character(ref_config) && ref_config == "auto")) {
  cat("\nReference: AUTO\n")
  ref_barcodes <- NULL
} else {
  ref_clusters <- as.character(ref_config)
  ref_barcodes <- colnames(sobj)[sobj$seurat_clusters %in% ref_clusters]
  cat(sprintf("\nReference: clusters %s → %d cells (%.1f%%)\n",
              paste(ref_clusters, collapse = ","),
              length(ref_barcodes),
              length(ref_barcodes) / ncol(sobj) * 100))

  if (length(ref_barcodes) < 100) {
    cat("  ⚠️  < 100 ref cells — consider adding clusters or 'auto'\n")
  }
}

# ============================================================================
# 6. Run SCEVAN
# ============================================================================

cat("\n══════════════════════════════════════════════════\n")
cat(sprintf("   Running SCEVAN pipelineCNA  [%s]\n", format(Sys.time(), "%H:%M:%S")))
cat("══════════════════════════════════════════════════\n\n")

t0 <- Sys.time()
old_wd <- getwd()
setwd(local_out)  # SCEVAN writes files to working dir

scevan_ok   <- FALSE
results     <- NULL
source_path <- "direct"  # will be "recovery" if crash path taken

tryCatch({
  if (is.null(ref_barcodes)) {
    results <- SCEVAN::pipelineCNA(
      count_mtx  = counts,
      sample     = scevan_cfg$sample_name,
      organism   = scevan_cfg$organism,
      par_cores  = scevan_cfg$par_cores,
      SUBCLONES  = isTRUE(scevan_cfg$subclones)
    )
  } else {
    results <- SCEVAN::pipelineCNA(
      count_mtx  = counts,
      sample     = scevan_cfg$sample_name,
      organism   = scevan_cfg$organism,
      par_cores  = scevan_cfg$par_cores,
      norm_cell  = ref_barcodes,
      SUBCLONES  = isTRUE(scevan_cfg$subclones)
    )
  }
  scevan_ok <- TRUE
  setwd(old_wd)

}, error = function(e) {
  setwd(old_wd)
  crash_msg <- sprintf("SCEVAN pipelineCNA error: %s", e$message)
  cat(sprintf("\n❌ %s\n", crash_msg))
  writeLines(crash_msg, file.path(local_out, "scevan_crash.log"))

  # Check for partial artifacts
  rdata_files <- list.files(local_out, pattern = "\\.RData$",
                            full.names = TRUE, recursive = TRUE)

  if (length(rdata_files) > 0) {
    cat(sprintf("\n📦 PARTIAL SUCCESS: %d .RData file(s) found → launching recovery\n",
                length(rdata_files)))
    source_path <<- "recovery"

    recovery_script <- file.path(module_root, "scripts", "03_recover_scevan_artifacts.R")
    if (!file.exists(recovery_script)) {
      stop("Recovery script not found: ", recovery_script)
    }

    Sys.setenv(LUNGMETA_ROOT = proj_root)
    ret <- system2(
      command = file.path(R.home("bin"), "Rscript"),
      args = c(shQuote(recovery_script), shQuote(normalizePath(cfg_path))),
      stdout = "", stderr = ""
    )

    if (ret == 0) {
      cat("\n✅ Recovery completed\n")
      scevan_ok <<- TRUE
    } else {
      stop(sprintf("Recovery failed (exit %d)", ret))
    }
  } else {
    stop("HARD FAILURE: No .RData artifacts. Must re-run from scratch.")
  }
})

elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)

if (!scevan_ok) {
  stop("SCEVAN did not complete successfully")
}

# ============================================================================
# 7. Standardize Output (THE KEY SECTION)
#    Both paths converge here to produce identical output format
# ============================================================================

cat("\n========== Standardizing output ==========\n")

if (source_path == "direct") {
  # ---- Path A: SCEVAN completed directly ----
  cat("Source: direct SCEVAN completion\n")

  # Parse the results data.frame
  if (is.data.frame(results)) {
    scevan_df <- results
    if (!"barcode" %in% colnames(scevan_df)) {
      scevan_df$barcode <- rownames(scevan_df)
    }
  } else {
    out_csv <- list.files(local_out, pattern = "\\.csv$",
                          full.names = TRUE, recursive = TRUE)
    if (length(out_csv) > 0) {
      scevan_df <- fread(out_csv[1])
    } else {
      saveRDS(results, file.path(local_out, "scevan_raw_unparsed.rds"))
      stop("Cannot parse SCEVAN direct output")
    }
  }

  # Find class column
  class_col <- intersect(c("class", "scevan_class", "tumor_normal"),
                         colnames(scevan_df))[1]
  if (is.na(class_col)) stop("No class column in SCEVAN output: ",
                              paste(colnames(scevan_df), collapse = ", "))

  # Build label vector for ALL cells
  bc_col <- if ("barcode" %in% colnames(scevan_df)) "barcode" else {
    scevan_df$barcode <- rownames(scevan_df); "barcode"
  }

  label_map <- setNames(scevan_df[[class_col]], scevan_df[[bc_col]])
  # Normalize to "tumor" / "non_tumor"
  label_map[label_map == "normal"] <- "non_tumor"

  all_labels <- rep("non_tumor", ncol(sobj))
  names(all_labels) <- colnames(sobj)
  matched_bc <- intersect(names(label_map), colnames(sobj))
  all_labels[matched_bc] <- label_map[matched_bc]

  sobj$scevan_label  <- all_labels[colnames(sobj)]
  sobj$scevan_source <- "direct"

  # Subclone if available
  subcl_col <- intersect(c("subclone", "Subclone"), colnames(scevan_df))
  if (length(subcl_col) > 0) {
    subcl_map <- setNames(scevan_df[[subcl_col[1]]], scevan_df[[bc_col]])
    sobj$scevan_subclone <- subcl_map[colnames(sobj)]
  } else {
    sobj$scevan_subclone <- NA_character_
  }

  # Save full SCEVAN result
  fwrite(scevan_df, file.path(local_out, "scevan_full_results.csv"))

} else {
  # ---- Path B: Recovery path ----
  cat("Source: recovery from .RData artifacts\n")

  # 03_recover already saved seurat_with_scevan_recovered.rds
  recovered_path <- file.path(main_out, "seurat_with_scevan_recovered.rds")
  if (!file.exists(recovered_path)) {
    recovered_path <- file.path(local_out, "seurat_with_scevan_recovered.rds")
  }
  if (!file.exists(recovered_path)) {
    stop("Recovery output not found: ", recovered_path)
  }

  sobj <- readRDS(recovered_path)

  # Verify required columns exist
  if (!"scevan_label" %in% colnames(sobj@meta.data)) {
    stop("scevan_label not in recovered object metadata")
  }

  sobj$scevan_source <- "recovery"

  if (!"scevan_subclone" %in% colnames(sobj@meta.data)) {
    sobj$scevan_subclone <- NA_character_
  }
}

# ============================================================================
# 8. Save Standardized Output (identical format regardless of path)
# ============================================================================

cat(sprintf("\nSaving standardized output → %s\n", basename(FINAL_OUTPUT)))
saveRDS(sobj, FINAL_OUTPUT)

# Quick summary
n_tumor    <- sum(sobj$scevan_label == "tumor", na.rm = TRUE)
n_nontumor <- sum(sobj$scevan_label != "tumor", na.rm = TRUE)

cat("\n")
cat("╔══════════════════════════════════════════════════════╗\n")
cat("║                01_scevan.R Complete                  ║\n")
cat(sprintf("║  Elapsed:  %.1f min                                  ║\n", elapsed))
cat(sprintf("║  Source:   %-42s ║\n", source_path))
cat(sprintf("║  Cells:    %6d                                    ║\n", ncol(sobj)))
cat(sprintf("║  Tumor:    %6d  (%.1f%%)                           ║\n",
            n_tumor, n_tumor / ncol(sobj) * 100))
cat(sprintf("║  Non-tumor:%6d  (%.1f%%)                           ║\n",
            n_nontumor, n_nontumor / ncol(sobj) * 100))
cat("╠══════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Output: %s\n", FINAL_OUTPUT))
cat("║  Next:   Rscript scripts/02_validate_scevan.R       ║\n")
cat("╚══════════════════════════════════════════════════════╝\n")
