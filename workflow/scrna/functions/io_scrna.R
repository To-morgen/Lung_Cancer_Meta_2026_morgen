# ============================================================================
# io_scrna.R — scRNA-seq I/O functions
# Layer: 2 (modality-level, scRNA-specific)
# Depends on: scripts/utils/utils_io.R (Layer 1)
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(Matrix)
})

# Load Layer 1
source(here("scripts", "utils", "utils_io.R"))

# ---- Cell Ranger path helpers ----

#' Standard Cell Ranger output paths for a sample
#' @param sample_id character
#' @param subdir character, default "full"
#' @return named list of paths
cellranger_paths <- function(sample_id, subdir = "full") {
  base <- file.path(get_ds_raw(), "cellranger_out", subdir, sample_id, "outs")
  list(
    base           = base,
    filtered_h5    = file.path(base, "filtered_feature_bc_matrix.h5"),
    raw_h5         = file.path(base, "raw_feature_bc_matrix.h5"),
    filtered_dir   = file.path(base, "filtered_feature_bc_matrix"),
    raw_dir        = file.path(base, "raw_feature_bc_matrix"),
    clusters       = file.path(base, "analysis", "clustering",
                               "gene_expression_graphclust", "clusters.csv"),
    metrics        = file.path(base, "metrics_summary.csv"),
    web_summary    = file.path(base, "web_summary.html")
  )
}

#' Standard scRNA results output paths
#' @param phase character, e.g. "02_qc"
#' @return named list of directories
scrna_output_dirs <- function(phase = "02_qc") {
  base <- here("results", "scrna", phase)
  dirs <- list(
    base        = base,
    soupx       = file.path(base, "soupx"),
    raw_seurat  = file.path(base, "raw_seurat"),
    qc_filtered = file.path(base, "qc_filtered"),
    doublets    = file.path(base, "doublets"),
    clean       = file.path(base, "clean"),
    plots       = file.path(base, "plots"),
    reports     = file.path(base, "reports")
  )
  # Ensure all dirs exist
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  return(dirs)
}

# ---- Matrix loaders ----

#' Load Cell Ranger filtered matrix (h5 or directory)
#' @param sample_id character
#' @param subdir character
#' @return dgCMatrix
load_filtered_matrix <- function(sample_id, subdir = "full") {
  paths <- cellranger_paths(sample_id, subdir)
  
  if (file.exists(paths$filtered_h5)) {
    log_msg(sprintf("Loading filtered h5: %s", sample_id))
    mat <- Seurat::Read10X_h5(paths$filtered_h5)
  } else if (dir.exists(paths$filtered_dir)) {
    log_msg(sprintf("Loading filtered dir: %s", sample_id))
    mat <- Seurat::Read10X(paths$filtered_dir)
  } else {
    stop(sprintf("No filtered matrix for %s at %s", sample_id, paths$base))
  }
  
  if (is.list(mat)) mat <- mat[["Gene Expression"]]
  log_msg(sprintf("%s: %d genes x %d cells", sample_id, nrow(mat), ncol(mat)))
  return(mat)
}

#' Load Cell Ranger raw (unfiltered) matrix
#' @param sample_id character
#' @param subdir character
#' @return dgCMatrix
load_raw_matrix <- function(sample_id, subdir = "full") {
  paths <- cellranger_paths(sample_id, subdir)
  
  if (file.exists(paths$raw_h5)) {
    log_msg(sprintf("Loading raw h5: %s", sample_id))
    mat <- Seurat::Read10X_h5(paths$raw_h5)
  } else if (dir.exists(paths$raw_dir)) {
    log_msg(sprintf("Loading raw dir: %s", sample_id))
    mat <- Seurat::Read10X(paths$raw_dir)
  } else {
    stop(sprintf("No raw matrix for %s at %s", sample_id, paths$base))
  }
  
  if (is.list(mat)) mat <- mat[["Gene Expression"]]
  log_msg(sprintf("%s raw: %d genes x %d barcodes", sample_id, nrow(mat), ncol(mat)))
  return(mat)
}

#' Load Cell Ranger clustering assignments
#' @param sample_id character
#' @return named vector or NULL
load_clusters <- function(sample_id, subdir = "full") {
  paths <- cellranger_paths(sample_id, subdir)
  
  if (!file.exists(paths$clusters)) {
    log_msg(sprintf("%s: no cluster file, SoupX will use automatic mode", sample_id), "warn")
    return(NULL)
  }
  
  cl <- data.table::fread(paths$clusters)
  clusters <- setNames(cl$Cluster, cl$Barcode)
  log_msg(sprintf("%s: %d cells in %d clusters", sample_id, length(clusters), length(unique(clusters))))
  return(clusters)
}

cat("[init] workflow/scrna/functions/io_scrna.R loaded\n")
