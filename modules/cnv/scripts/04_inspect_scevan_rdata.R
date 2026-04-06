#!/usr/bin/env Rscript
# DEBUG/INSPECTION ONLY — not called by pipeline

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(data.table)
  library(Seurat)
})

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

load_rdata_as_list <- function(path) {
  e <- new.env(parent = emptyenv())
  nm <- load(path, envir = e)
  stats::setNames(lapply(nm, function(x) get(x, envir = e)), nm)
}

safe_dim <- function(x) {
  d <- dim(x)
  if (is.null(d)) c(NA_integer_, NA_integer_) else d
}

has_cols <- function(x, candidates) {
  if (!is.data.frame(x)) return(FALSE)
  any(candidates %in% colnames(x))
}

count_overlap <- function(v, target) {
  if (is.null(v)) return(0L)
  sum(as.character(v) %in% target)
}

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

module_root <- here::here()
proj_root <- Sys.getenv("LUNGMETA_ROOT")
if (proj_root == "") {
  proj_root <- normalizePath(file.path(module_root, "..", ".."))
}

cfg <- yaml::read_yaml(file.path(module_root, "configs", "cnv_params.yaml"))

artifact_dir <- file.path(module_root, cfg$output$module_results, "scevan", "output")
if (!dir.exists(artifact_dir)) {
  artifact_dir <- file.path(module_root, cfg$output$module_results, "scevan")
}

input_path <- file.path(proj_root, cfg$input$seurat_object)

cat("Artifact dir: ", artifact_dir, "\n", sep = "")
cat("Input Seurat: ", input_path, "\n", sep = "")

# -----------------------------------------------------------------------------
# Load mother object and reference cells
# -----------------------------------------------------------------------------

sobj <- readRDS(input_path)
cell_barcodes <- colnames(sobj)

cluster_col <- if (!is.null(cfg$input$cluster_column)) cfg$input$cluster_column else "seurat_clusters"
ref_clusters <- as.character(cfg$scevan$reference_clusters)
ref_barcodes <- colnames(sobj)[as.character(sobj[[cluster_col]][, 1]) %in% ref_clusters]

cat("Total cells: ", length(cell_barcodes), "\n", sep = "")
cat("Reference cells: ", length(ref_barcodes), "\n", sep = "")
cat("Reference clusters: ", paste(ref_clusters, collapse = ", "), "\n", sep = "")

# -----------------------------------------------------------------------------
# Inspect all RData objects
# -----------------------------------------------------------------------------

rdata_files <- list.files(artifact_dir, pattern = "\\.RData$", full.names = TRUE)
if (length(rdata_files) == 0) stop("No .RData files found in: ", artifact_dir)

out <- list()
k <- 1L

for (f in rdata_files) {
  obj_list <- load_rdata_as_list(f)

  for (nm in names(obj_list)) {
    obj <- obj_list[[nm]]
    d <- safe_dim(obj)

    row_ov_all <- count_overlap(rownames(obj), cell_barcodes)
    col_ov_all <- count_overlap(colnames(obj), cell_barcodes)
    row_ov_ref <- count_overlap(rownames(obj), ref_barcodes)
    col_ov_ref <- count_overlap(colnames(obj), ref_barcodes)

    bc_col_present <- if (is.data.frame(obj)) {
      paste(intersect(
        c("barcode", "Barcode", "cell", "Cell", "cell_id", "CellID", "X"),
        colnames(obj)
      ), collapse = ";")
    } else {
      ""
    }

    label_col_present <- if (is.data.frame(obj)) {
      paste(intersect(
        c("class", "scevan_class", "tumor_normal", "prediction", "status"),
        colnames(obj)
      ), collapse = ";")
    } else {
      ""
    }

    subclone_col_present <- if (is.data.frame(obj)) {
      paste(intersect(
        c("subclone", "Subclone", "SUBCLONE"),
        colnames(obj)
      ), collapse = ";")
    } else {
      ""
    }

    out[[k]] <- data.frame(
      file = basename(f),
      object = nm,
      class = paste(class(obj), collapse = ";"),
      nrow = d[1],
      ncol = d[2],
      row_overlap_all = row_ov_all,
      col_overlap_all = col_ov_all,
      row_overlap_ref = row_ov_ref,
      col_overlap_ref = col_ov_ref,
      barcode_columns = bc_col_present,
      label_columns = label_col_present,
      subclone_columns = subclone_col_present,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}

manifest <- rbindlist(out, fill = TRUE)

manifest[, max_cell_overlap := pmax(row_overlap_all, col_overlap_all, na.rm = TRUE)]
manifest[, max_ref_overlap := pmax(row_overlap_ref, col_overlap_ref, na.rm = TRUE)]
manifest[, has_label_cols := label_columns != ""]
manifest[, has_subclone_cols := subclone_columns != ""]

setorder(manifest, -has_label_cols, -has_subclone_cols, -max_cell_overlap, -max_ref_overlap)

out_csv <- file.path(proj_root, cfg$output$main_results, "reports", "scevan_rdata_deep_inspection.csv")
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
fwrite(manifest, out_csv)

cat("\nTop candidates:\n")
print(head(manifest, 30))

cat("\nSaved to:\n", out_csv, "\n", sep = "")
