# ============================================================================
# io_scrna.R — I/O utilities for scRNA-seq pipeline
#
# Provides:
#   - load_sample_list()      → character vector of sample IDs
#   - load_sample_groups()    → named vector: sample → group
#   - scrna_output_dirs()     → phase-specific directory list
#   - load_dataset_config()   → full dataset YAML
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(yaml)
})

source(here("scripts", "utils", "utils_io.R"))

# ---- Dataset Config ----

load_dataset_config <- function() {
  cfg_path <- Sys.getenv("DS_CONFIG")
  if (cfg_path == "" || !file.exists(cfg_path)) {
    cfg_path <- here("configs", "datasets",
                     paste0(Sys.getenv("DS_NAME", "UNKNOWN"), ".yaml"))
  }
  if (!file.exists(cfg_path)) stop("Dataset config not found: ", cfg_path)
  cfg <- yaml::read_yaml(cfg_path)
  log_msg(sprintf("Loaded dataset config: %s", cfg$dataset_id))
  cfg
}

# ---- Sample List & Groups ----

load_sample_list <- function() {
  cfg <- load_dataset_config()
  sids <- names(cfg$samples)
  log_msg(sprintf("Samples: %s", paste(sids, collapse = ", ")))
  sids
}

load_sample_groups <- function() {
  cfg <- load_dataset_config()
  grp <- sapply(cfg$samples, function(s) s$group)
  names(grp) <- names(cfg$samples)
  grp
}

# ---- Output Directories (phase-aware) ----

scrna_output_dirs <- function(phase = "02_qc") {
  base <- here("results", "scrna", phase)

  if (phase == "02_qc") {
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
  } else if (phase == "03_normalize") {
    dirs <- list(
      base       = base,
      merged     = file.path(base, "merged"),
      normalized = file.path(base, "normalized"),
      pca        = file.path(base, "pca"),
      plots      = file.path(base, "plots"),
      reports    = file.path(base, "reports")
    )
  } else if (phase == "04_integrate") {
    dirs <- list(
      base     = base,
      harmony  = file.path(base, "harmony"),
      clusters = file.path(base, "clusters"),
      umap     = file.path(base, "umap"),
      plots    = file.path(base, "plots"),
      reports  = file.path(base, "reports")
    )
  } else {
    # Generic fallback
    dirs <- list(
      base    = base,
      plots   = file.path(base, "plots"),
      reports = file.path(base, "reports")
    )
  }

  # Ensure all dirs exist
  lapply(dirs, function(d) dir.create(d, recursive = TRUE, showWarnings = FALSE))
  invisible(dirs)
}
