# ============================================================================
# qc_utils.R — scRNA-seq QC utility functions
# Layer: 2 (modality-level, scRNA-specific)
# Depends on: scripts/utils/utils_io.R (Layer 1)
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(yaml)
})

source(here("scripts", "utils", "utils_io.R"))

# ====================================================================
# QC Parameters (YAML config chain)
# ====================================================================

#' Load QC parameters with priority chain:
#'   dataset private config → project params → hardcoded defaults
load_qc_params <- function() {
  # Hardcoded defaults (safety net)
  defaults <- list(
    species = "mouse",
    filter_method = "MAD",
    n_mad = 3,
    min_genes = 200,
    max_genes = Inf,
    min_umi = 500,
    max_mt_pct = 20,
    min_cells_per_gene = 3,
    mito_pattern = "^mt-",
    ribo_pattern = "^Rp[sl]",
    doublet_seed = 42
  )

  # Project-level config
  proj_path <- here("configs", "params", "scrna_qc_params.yaml")
  if (file.exists(proj_path)) {
    proj <- yaml::read_yaml(proj_path)
    if (!is.null(proj$cell_filter)) {
      for (k in names(proj$cell_filter)) {
        defaults[[k]] <- proj$cell_filter[[k]]
      }
    }
    # Read top-level species from project config
    if (!is.null(proj$species)) defaults$species <- proj$species
    # Pass through additional top-level sections
    if (!is.null(proj$normalization)) defaults$normalization <- proj$normalization
    if (!is.null(proj$cell_cycle)) defaults$cell_cycle <- proj$cell_cycle
    if (!is.null(proj$clustering))    defaults$clustering    <- proj$clustering
    if (!is.null(proj$integration))   defaults$integration   <- proj$integration
    log_msg("QC params: loaded project config")
  }

  # Dataset-level overrides
  tryCatch({
    ds_cfg <- load_dataset_config()
    if (!is.null(ds_cfg$qc_overrides)) {
      for (k in names(ds_cfg$qc_overrides)) {
        defaults[[k]] <- ds_cfg$qc_overrides[[k]]
      }
      log_msg("QC params: dataset overrides applied")
    }
    # Top-level species from dataset config (highest priority)
    if (!is.null(ds_cfg$species)) {
      defaults$species <- ds_cfg$species
      log_msg(sprintf("QC params: species = '%s' (from dataset config)", ds_cfg$species))
    }
  }, error = function(e) {
    log_msg(sprintf("QC params: no dataset overrides (%s)", e$message), "WARN")
  })

  # Derive patterns from species
  if (defaults$species == "mouse") {
    defaults$mito_pattern <- "^mt-"
    defaults$ribo_pattern <- "^Rp[sl]"
  } else if (defaults$species == "human") {
    defaults$mito_pattern <- "^MT-"
    defaults$ribo_pattern <- "^RP[SL]"
  }

  log_msg(sprintf("  species=%s, method=%s, n_mad=%d, hard_min_genes=%d, hard_max_mt=%d%%, mito_pattern=%s",
                   defaults$species, defaults$filter_method, defaults$n_mad,
                   defaults$min_genes, defaults$max_mt_pct, defaults$mito_pattern))
  defaults
}


#' Load sample list from dataset config
load_sample_list <- function() {
  cfg <- load_dataset_config()
  ids <- names(cfg$samples)
  log_msg(sprintf("Samples: %s", paste(ids, collapse = ", ")))
  ids
}

#' Get sample-to-group mapping
load_sample_groups <- function() {
  cfg <- load_dataset_config()
  sapply(cfg$samples, function(s) s$group)
}

# ====================================================================
# QC Metrics
# ====================================================================

add_mito_pct <- function(sobj, pattern = "^mt-") {
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = pattern)
  sobj
}

add_ribo_pct <- function(sobj, pattern = "^Rp[sl]") {
  sobj[["percent.ribo"]] <- PercentageFeatureSet(sobj, pattern = pattern)
  sobj
}

add_hemo_pct <- function(sobj, species = "mouse") {
  pattern <- switch(species, mouse = "^Hb[ab]-", human = "^HB[AB]", "^Hb")
  sobj[["percent.hb"]] <- PercentageFeatureSet(sobj, pattern = pattern)
  sobj
}

# ====================================================================
# MAD-based Outlier Detection
# ====================================================================

#' Detect outliers using MAD method
#'
#' @param values numeric vector
#' @param n_mad number of MADs (default 3)
#' @param log_transform apply log10 before computing (for right-skewed data)
#' @param direction "both" | "upper" | "lower"
#' @param hard_min absolute minimum (in original space)
#' @param hard_max absolute maximum (in original space)
#' @return logical vector (TRUE = outlier = should REMOVE)
detect_outliers_mad <- function(values,
                                 n_mad         = 3,
                                 log_transform = FALSE,
                                 direction     = "both",
                                 hard_min      = -Inf,
                                 hard_max      = Inf) {

  x <- if (log_transform) log10(values + 1) else values

  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, na.rm = TRUE)

  # Compute bounds in (possibly log) space
  lower_x <- if (direction %in% c("lower", "both")) med - n_mad * mad_val else -Inf
  upper_x <- if (direction %in% c("upper", "both")) med + n_mad * mad_val else Inf

  # Transform back if log
  if (log_transform) {
    lower_orig <- 10^lower_x - 1
    upper_orig <- 10^upper_x - 1
  } else {
    lower_orig <- lower_x
    upper_orig <- upper_x
  }

  # Apply hard bounds
  lower_final <- max(lower_orig, hard_min)
  upper_final <- min(upper_orig, hard_max)

  is_outlier <- (values < lower_final) | (values > upper_final)

  # Return with attributes for reporting
  attr(is_outlier, "lower") <- lower_final
  attr(is_outlier, "upper") <- upper_final
  attr(is_outlier, "median") <- if (log_transform) 10^med - 1 else med
  attr(is_outlier, "mad") <- mad_val
  attr(is_outlier, "n_outlier") <- sum(is_outlier, na.rm = TRUE)

  is_outlier
}

#' Apply full MAD-based QC to a Seurat object
#'
#' @param sobj Seurat object (with percent.mt already added)
#' @param params list from load_qc_params()
#' @return list(pass = logical, report = data.frame)
mad_qc_filter <- function(sobj, params) {

  n_mad <- params$n_mad

  # nFeature_RNA — log10 space, both directions
  out_nf <- detect_outliers_mad(
    sobj$nFeature_RNA,
    n_mad = n_mad, log_transform = TRUE, direction = "both",
    hard_min = params$min_genes
  )

  # nCount_RNA — log10 space, both directions
  out_nc <- detect_outliers_mad(
    sobj$nCount_RNA,
    n_mad = n_mad, log_transform = TRUE, direction = "both",
    hard_min = params$min_umi
  )

  # percent.mt — linear space, upper only
  out_mt <- detect_outliers_mad(
    sobj$percent.mt,
    n_mad = n_mad, log_transform = FALSE, direction = "upper",
    hard_max = params$max_mt_pct
  )

  # Combined
  is_outlier <- out_nf | out_nc | out_mt
  pass <- !is_outlier

  # Report
  report <- data.frame(
    metric       = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    space        = c("log10", "log10", "linear"),
    direction    = c("both", "both", "upper"),
    lower_thresh = c(attr(out_nf, "lower"), attr(out_nc, "lower"), NA),
    upper_thresh = c(attr(out_nf, "upper"), attr(out_nc, "upper"), attr(out_mt, "upper")),
    median_val   = c(attr(out_nf, "median"), attr(out_nc, "median"), attr(out_mt, "median")),
    n_outlier    = c(attr(out_nf, "n_outlier"), attr(out_nc, "n_outlier"), attr(out_mt, "n_outlier")),
    stringsAsFactors = FALSE
  )

  list(pass = pass, report = report)
}

# ====================================================================
# Summary helpers
# ====================================================================

qc_summary_table <- function(sobj, sample_id = "sample") {
  md <- sobj@meta.data
  data.frame(
    sample          = sample_id,
    n_cells         = nrow(md),
    median_genes    = median(md$nFeature_RNA),
    median_umi      = median(md$nCount_RNA),
    median_mt_pct   = round(median(md$percent.mt, na.rm = TRUE), 2),
    mean_mt_pct     = round(mean(md$percent.mt, na.rm = TRUE), 2),
    stringsAsFactors = FALSE
  )
}

cat("[init] workflow/scrna/functions/qc_utils.R loaded\n")
