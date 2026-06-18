# ============================================================================
# utils_io.R — Project-level I/O and path utilities
# Layer: 1 (project-wide, no modality-specific logic)
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(yaml)
})

resolve_project_root <- function(start = getwd()) {
  env_root <- Sys.getenv("PROJECT_ROOT", unset = Sys.getenv("LUNGMETA_ROOT", unset = ""))
  if (nzchar(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = FALSE))
  }

  candidates <- unique(c(
    normalizePath(start, winslash = "/", mustWork = FALSE),
    normalizePath(here::here(), winslash = "/", mustWork = FALSE)
  ))

  is_project_root <- function(path) {
    file.exists(file.path(path, "Lung_Cancer_Meta_2026.Rproj")) &&
      dir.exists(file.path(path, "scripts")) &&
      dir.exists(file.path(path, "configs"))
  }

  ascend_until_root <- function(path) {
    current <- path
    if (file.exists(current) && !dir.exists(current)) {
      current <- dirname(current)
    }
    current <- normalizePath(current, winslash = "/", mustWork = FALSE)

    repeat {
      if (is_project_root(current)) return(current)
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
    NULL
  }

  for (candidate in candidates) {
    root <- ascend_until_root(candidate)
    if (!is.null(root)) return(root)
  }

  stop("Could not resolve project root. Set PROJECT_ROOT or run from inside the repository.")
}

project_path <- function(...) {
  file.path(resolve_project_root(), ...)
}

#' Load .env.sh variables into R session
load_env <- function(env_file = project_path(".env.sh")) {
  if (!file.exists(env_file)) stop(".env.sh not found: ", env_file)
  lines <- readLines(env_file)
  lines <- lines[grepl("^export ", lines)]
  for (line in lines) {
    line <- sub("^export ", "", line)
    eq_pos <- regexpr("=", line, fixed = TRUE)
    key <- substr(line, 1, eq_pos - 1)
    val <- substr(line, eq_pos + 1, nchar(line))
    val <- gsub('^"|"$', '', val)
    # Expand variables
    val <- gsub("\\$HOME", Sys.getenv("HOME"), val)
    val <- gsub("\\$\\{HOME\\}", Sys.getenv("HOME"), val)
    # Expand ${VAR} references to already-set vars
    while (grepl("\\$\\{[A-Z_]+\\}", val)) {
      m <- regmatches(val, regexpr("\\$\\{[A-Z_]+\\}", val))
      var_name <- gsub("\\$\\{|\\}", "", m)
      val <- sub(m, Sys.getenv(var_name, ""), val, fixed = TRUE)
    }
    do.call(Sys.setenv, setNames(list(val), key))
  }
  invisible(NULL)
}

#' Simple timestamped logger
#' @param msg character message
#' @param level "info", "warn", "error" (default "info")
log_msg <- function(msg, level = "info") {
  timestamp <- format(Sys.time(), "[%H:%M:%S]")
  prefix <- switch(tolower(level),
    warn  = "[warn]",
    error = "[ERROR]",
    "[info]"
  )
  cat(sprintf("%s %s %s\n", timestamp, prefix, msg))
}

#' Load dataset-specific config
#' Priority: DS_CONFIG env var (Snakemake) > DATASET_ID env var (legacy)
#' @param dataset_id character (optional, auto-detect from env)
#' @return list (parsed YAML)
load_dataset_config <- function(dataset_id = NULL) {
  # ── Priority 1: DS_CONFIG env var (new standard, used by Snakemake) ──
  cfg_path <- Sys.getenv("DS_CONFIG", unset = "")

  if (nzchar(cfg_path)) {
    # Resolve relative paths via here()
    if (!startsWith(cfg_path, "/")) {
      cfg_path <- here(cfg_path)
    }
    if (file.exists(cfg_path)) {
      cfg <- yaml::read_yaml(cfg_path)
      log_msg(sprintf("Loaded dataset config: %s (from %s)", cfg$dataset_id, cfg_path))
      return(cfg)
    } else {
      log_msg(sprintf("DS_CONFIG set but file not found: %s", cfg_path), "WARN")
    }
  }

  # ── Priority 2: DATASET_ID env var (legacy, backward-compatible) ──
  if (is.null(dataset_id)) {
    dataset_id <- Sys.getenv("DATASET_ID", unset = "")
    if (!nzchar(dataset_id)) {
      load_env()
      dataset_id <- Sys.getenv("DATASET_ID", unset = "")
    }
  }

  if (!nzchar(dataset_id)) stop("Neither DS_CONFIG nor DATASET_ID is set")

  # Legacy path: configs_private/datasets/{DATASET_ID}.yaml
  cfg_path <- project_path("configs_private", "datasets", paste0(dataset_id, ".yaml"))
  if (!file.exists(cfg_path)) {
    # Also try configs/datasets/
    cfg_path <- project_path("configs", "datasets", paste0(dataset_id, ".yaml"))
  }
  if (!file.exists(cfg_path)) stop("Config not found for dataset_id: ", dataset_id)

  cfg <- yaml::read_yaml(cfg_path)
  log_msg(sprintf("Loaded dataset config: %s (from %s)", cfg$dataset_id, cfg_path))
  cfg
}

cat("[init] scripts/utils/utils_io.R loaded\n")
