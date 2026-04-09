# ============================================================================
# io_scrna.R — I/O utilities for scRNA-seq pipeline
#
# Provides:
#   - load_dataset_config()   → full dataset YAML
#   - load_sample_list()      → character vector of sample IDs
#   - load_sample_groups()    → named vector: sample → group
#   - scrna_output_dirs()     → phase-specific directory list (generic)
#   - scrna_integrate_dirs()  → Phase 04 dirs
#   - scrna_cluster_dirs()    → Phase 05 dirs
#   - scrna_annotate_dirs()   → Phase 06 dirs  [NEW]
#   - scrna_de_dirs()         → Phase 08 dirs  [NEW]
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(yaml)
})

source(here("scripts", "utils", "utils_io.R"))

# ---- Internal helper: resolve base path with DS_PREFIX ----
.scrna_base <- function(phase) {
  ds_prefix <- Sys.getenv("DS_PREFIX", unset = "")
  if (ds_prefix != "") {
    base <- here("results", "scrna", ds_prefix, phase)
  } else {
    base <- here("results", "scrna", phase)
  }
  base
}

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

# ---- Output Directories (phase-aware, generic) ----

scrna_output_dirs <- function(phase = "02_qc") {
  base <- .scrna_base(phase)

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
  } else if (phase == "05_cluster") {
    dirs <- list(
      base    = base,
      objects = file.path(base, "objects"),
      plots   = file.path(base, "plots"),
      reports = file.path(base, "reports")
    )
  } else if (phase == "06_annotate") {
    dirs <- list(
      base    = base,
      markers = file.path(base, "markers"),
      singler = file.path(base, "singler"),
      objects = file.path(base, "objects"),
      plots   = file.path(base, "plots"),
      reports = file.path(base, "reports")
    )
  } else if (phase == "08_de") {
    dirs <- list(
      base       = base,
      pseudobulk = file.path(base, "pseudobulk"),
      deseq2     = file.path(base, "deseq2"),
      enrichment = file.path(base, "enrichment"),
      plots      = file.path(base, "plots"),
      reports    = file.path(base, "reports")
    )
  } else {
    # Generic fallback
    dirs <- list(
      base    = base,
      plots   = file.path(base, "plots"),
      reports = file.path(base, "reports")
    )
  }

  lapply(dirs, function(d) dir.create(d, recursive = TRUE, showWarnings = FALSE))
  ds_prefix <- Sys.getenv("DS_PREFIX", unset = "")
  if (ds_prefix != "") {
    log_msg(sprintf("Output dirs: results/scrna/%s/%s (%d subdirs)",
                    ds_prefix, phase, length(dirs)))
  }
  invisible(dirs)
}

# ---- Phase-specific convenience helpers ----

scrna_integrate_dirs <- function() {
  base <- .scrna_base("04_integrate")
  dirs <- list(
    base    = base,
    harmony = file.path(base, "harmony"),
    plots   = file.path(base, "plots"),
    reports = file.path(base, "reports")
  )
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dirs
}

scrna_cluster_dirs <- function() {
  base <- .scrna_base("05_cluster")
  dirs <- list(
    base    = base,
    objects = file.path(base, "objects"),
    plots   = file.path(base, "plots"),
    reports = file.path(base, "reports")
  )
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dirs
}

scrna_annotate_dirs <- function() {
  base <- .scrna_base("06_annotate")
  dirs <- list(
    base    = base,
    markers = file.path(base, "markers"),
    singler = file.path(base, "singler"),
    objects = file.path(base, "objects"),
    plots   = file.path(base, "plots"),
    reports = file.path(base, "reports")
  )
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dirs
}

scrna_de_dirs <- function() {
  base <- .scrna_base("08_de")
  dirs <- list(
    base       = base,
    pseudobulk = file.path(base, "pseudobulk"),
    deseq2     = file.path(base, "deseq2"),
    enrichment = file.path(base, "enrichment"),
    plots      = file.path(base, "plots"),
    reports    = file.path(base, "reports")
  )
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dirs
}

cat("[init] io_scrna.R loaded (DS_PREFIX-aware, phases 02–08)\n")
