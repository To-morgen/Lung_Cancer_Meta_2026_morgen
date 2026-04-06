#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(here)
  library(yaml)
})

# =============================================================================
# 03_recover_scevan_artifacts.R
#
# Purpose:
#   Safely recover tumor membership from native SCEVAN artifacts after a
#   late-stage crash.
#
# Recovery policy:
#   - Promote all raw artifacts to project-level result directory
#   - Recover tumor membership ONLY from results.com
#   - Reject biologically implausible recovery automatically
#
# Biological note:
#   results.com behaves like a tumor-focused subclone object in this run.
#   CNA_mtx_relat is NOT treated as a tumor-membership carrier.
# =============================================================================

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

#' Read config with lightweight fallback.
#'
#' Args:
#'   cfg_path: YAML file path
#'
#' Returns:
#'   Parsed config list
read_cfg <- function(cfg_path) {
  if (!file.exists(cfg_path)) stop("Config not found: ", cfg_path)
  yaml::read_yaml(cfg_path)
}

#' Get a nested config value using multiple candidate paths.
#'
#' Args:
#'   cfg: Parsed config list
#'   paths: List of character vectors, each a candidate nested path
#'   default: Default value if all paths fail
#'
#' Returns:
#'   Retrieved value or default
get_cfg_value <- function(cfg, paths, default = NULL) {
  for (p in paths) {
    x <- cfg
    ok <- TRUE
    for (nm in p) {
      if (!is.list(x) || is.null(x[[nm]])) {
        ok <- FALSE
        break
      }
      x <- x[[nm]]
    }
    if (ok) return(x)
  }
  default
}

#' Load all objects from an .RData file.
#'
#' Args:
#'   path: .RData path
#'
#' Returns:
#'   Named list of objects
load_rdata_as_list <- function(path) {
  e <- new.env(parent = emptyenv())
  nm <- load(path, envir = e)
  stats::setNames(lapply(nm, function(x) get(x, envir = e)), nm)
}

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

module_root <- here::here()
proj_root <- Sys.getenv("LUNGMETA_ROOT")
if (proj_root == "") {
  proj_root <- normalizePath(file.path(module_root, "..", ".."))
}

cfg <- read_cfg(file.path(module_root, "configs", "cnv_params.yaml"))

input_path <- file.path(
  proj_root,
  get_cfg_value(cfg, list(c("input", "seurat_object")), "results/scrna/06_annotate/objects/seurat_annotated.rds")
)

artifact_dir <- file.path(
  module_root,
  get_cfg_value(cfg, list(c("output", "module_results")), "results"),
  "scevan", "output"
)
if (!dir.exists(artifact_dir)) {
  artifact_dir <- file.path(
    module_root,
    get_cfg_value(cfg, list(c("output", "module_results")), "results"),
    "scevan"
  )
}

main_out <- file.path(
  proj_root,
  get_cfg_value(cfg, list(c("output", "main_results")), "results/scrna/10_cnv"),
  "scevan"
)
plot_out <- file.path(
  proj_root,
  get_cfg_value(cfg, list(c("output", "main_results")), "results/scrna/10_cnv"),
  "plots"
)
report_out <- file.path(
  proj_root,
  get_cfg_value(cfg, list(c("output", "main_results")), "results/scrna/10_cnv"),
  "reports"
)

dir.create(main_out, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_out, recursive = TRUE, showWarnings = FALSE)
dir.create(report_out, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Load mother object and curated references
# -----------------------------------------------------------------------------

sobj <- readRDS(input_path)
cell_barcodes <- colnames(sobj)

cluster_col <- get_cfg_value(
  cfg,
  list(c("input", "cluster_column"), c("input", "cluster_col")),
  "seurat_clusters"
)

reference_clusters <- as.character(get_cfg_value(
  cfg,
  list(
    c("input", "reference_clusters"),
    c("input", "ref_clusters"),
    c("scevan", "reference_clusters"),
    c("reference", "clusters")
  ),
  c("3", "4", "5", "6", "11", "14", "15")
))

if (!cluster_col %in% colnames(sobj@meta.data)) {
  stop("Cluster column not found in Seurat metadata: ", cluster_col)
}

cluster_vec <- as.character(sobj[[cluster_col]][, 1])
ref_barcodes <- colnames(sobj)[cluster_vec %in% reference_clusters]

cat("Project root: ", proj_root, "\n", sep = "")
cat("Artifact dir: ", artifact_dir, "\n", sep = "")
cat("Input Seurat: ", input_path, "\n", sep = "")
cat("Total cells: ", length(cell_barcodes), "\n", sep = "")
cat("Reference clusters: ", paste(reference_clusters, collapse = ", "), "\n", sep = "")
cat("Reference cells: ", length(ref_barcodes), "\n", sep = "")

# -----------------------------------------------------------------------------
# Promote raw artifacts
# -----------------------------------------------------------------------------

artifact_files <- list.files(artifact_dir, full.names = TRUE, recursive = FALSE)
copy_ok <- file.copy(artifact_files, main_out, overwrite = TRUE)

artifact_manifest <- data.frame(
  file = basename(artifact_files),
  src = artifact_files,
  dst = file.path(main_out, basename(artifact_files)),
  copied = copy_ok,
  stringsAsFactors = FALSE
)
fwrite(artifact_manifest, file.path(main_out, "scevan_artifact_manifest.csv"))

plot_files <- grep("\\.(png|pdf)$", artifact_files, value = TRUE, ignore.case = TRUE)
if (length(plot_files) > 0) {
  file.copy(plot_files, plot_out, overwrite = TRUE)
}

# -----------------------------------------------------------------------------
# Deep inspect RData files and save manifest
# -----------------------------------------------------------------------------

rdata_files <- list.files(artifact_dir, pattern = "\\.RData$", full.names = TRUE)
sig_list <- list()
k <- 1L

for (f in rdata_files) {
  obj_list <- load_rdata_as_list(f)
  for (nm in names(obj_list)) {
    obj <- obj_list[[nm]]
    nr <- if (!is.null(dim(obj))) dim(obj)[1] else NA_integer_
    nc <- if (!is.null(dim(obj))) dim(obj)[2] else NA_integer_
    sig_list[[k]] <- data.frame(
      file = basename(f),
      object = nm,
      class = paste(class(obj), collapse = ";"),
      nrow = nr,
      ncol = nc,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}

if (length(sig_list) > 0) {
  fwrite(rbindlist(sig_list), file.path(report_out, "scevan_rdata_object_manifest.csv"))
}

# -----------------------------------------------------------------------------
# Recover tumor membership ONLY from results.com
# -----------------------------------------------------------------------------

subclone_rdata <- file.path(artifact_dir, "LLC_tumor_CNAmtxSubclones.RData")
if (!file.exists(subclone_rdata)) {
  stop("Expected subclone RData not found: ", subclone_rdata)
}

obj_list <- load_rdata_as_list(subclone_rdata)
if (!"results.com" %in% names(obj_list)) {
  stop("results.com not found in: ", subclone_rdata)
}

results_com <- obj_list[["results.com"]]

if (is.null(colnames(results_com))) {
  stop("results.com has no colnames; cannot map tumor candidates.")
}

tumor_barcodes <- unique(colnames(results_com))
tumor_barcodes <- tumor_barcodes[tumor_barcodes %in% cell_barcodes]

candidate_n <- length(tumor_barcodes)
candidate_frac <- candidate_n / length(cell_barcodes)

ref_overlap_n <- sum(tumor_barcodes %in% ref_barcodes)
ref_overlap_frac <- if (length(ref_barcodes) > 0) ref_overlap_n / length(ref_barcodes) else NA_real_

cat("Tumor candidates from results.com: ", candidate_n, "\n", sep = "")
cat("Tumor candidate fraction: ", round(100 * candidate_frac, 1), "%\n", sep = "")
cat("Reference contamination: ", ref_overlap_n, " / ", length(ref_barcodes),
    " = ", round(100 * ref_overlap_frac, 3), "%\n", sep = "")

# Plausibility guards
if (candidate_n < 1000) {
  stop("Recovered tumor candidate set is too small; refusing promotion.")
}
if (candidate_frac > 0.9) {
  stop("Recovered tumor candidate fraction is implausibly high (>90%); refusing promotion.")
}
if (!is.na(ref_overlap_frac) && ref_overlap_frac > 0.05) {
  stop("Reference contamination exceeds 5%; refusing promotion.")
}

# -----------------------------------------------------------------------------
# Write recovered labels
# -----------------------------------------------------------------------------

full_labels <- data.frame(
  barcode = cell_barcodes,
  scevan_call = ifelse(cell_barcodes %in% tumor_barcodes, "tumor", "non_tumor"),
  recovery_mode = "tumor_membership_from_results.com",
  source_file = "LLC_tumor_CNAmtxSubclones.RData",
  source_object = "results.com",
  subclone = NA_character_,
  stringsAsFactors = FALSE
)

raw_labels <- full_labels[full_labels$scevan_call == "tumor", , drop = FALSE]

fwrite(raw_labels, file.path(main_out, "scevan_recovered_labels_raw.csv"))
fwrite(full_labels, file.path(main_out, "scevan_recovered_labels_full.csv"))
fwrite(full_labels, file.path(report_out, "scevan_recovered_labels_full.csv"))

meta_df <- full_labels[, c("barcode", "scevan_call", "recovery_mode"), drop = FALSE]
rownames(meta_df) <- meta_df$barcode
meta_df$barcode <- NULL
colnames(meta_df) <- c("scevan_call", "scevan_recovery_mode")

sobj_labeled <- AddMetaData(sobj, metadata = meta_df)
saveRDS(sobj_labeled, file.path(main_out, "seurat_with_scevan_recovered.rds"))

# -----------------------------------------------------------------------------
# Summaries
# -----------------------------------------------------------------------------

summary_dt <- data.table(
  barcode = cell_barcodes,
  cluster = cluster_vec,
  scevan_call = full_labels$scevan_call
)

cluster_summary <- summary_dt[, .(
  n_total = .N,
  n_tumor = sum(scevan_call == "tumor"),
  n_non_tumor = sum(scevan_call == "non_tumor"),
  pct_tumor = round(100 * sum(scevan_call == "tumor") / .N, 1)
), by = cluster][order(suppressWarnings(as.integer(cluster)), cluster)]

ref_check <- summary_dt[, .(
  n_total = .N,
  n_tumor = sum(scevan_call == "tumor"),
  pct_tumor = round(100 * sum(scevan_call == "tumor") / .N, 3)
), by = .(is_reference_cluster = cluster %in% reference_clusters)]

fwrite(cluster_summary, file.path(report_out, "scevan_recovered_tumor_by_cluster.csv"))
fwrite(ref_check, file.path(report_out, "scevan_recovered_reference_check.csv"))

run_note <- data.frame(
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  recovery_status = "success",
  candidate_n = candidate_n,
  candidate_frac = round(candidate_frac, 4),
  ref_overlap_n = ref_overlap_n,
  ref_overlap_frac = round(ref_overlap_frac, 6),
  source_file = "LLC_tumor_CNAmtxSubclones.RData",
  source_object = "results.com",
  stringsAsFactors = FALSE
)
fwrite(run_note, file.path(report_out, "scevan_recovery_run_note.csv"))

cat("\nRecovery succeeded.\n")
cat("Recovered labels file: ", file.path(main_out, "scevan_recovered_labels_full.csv"), "\n", sep = "")
cat("Recovered Seurat file: ", file.path(main_out, "seurat_with_scevan_recovered.rds"), "\n", sep = "")
cat("Reference check file: ", file.path(report_out, "scevan_recovered_reference_check.csv"), "\n", sep = "")
cat("Done.\n")
