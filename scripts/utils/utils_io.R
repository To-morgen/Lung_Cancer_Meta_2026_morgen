# ============================================================================
# utils_io.R — Project-level I/O and path utilities
# Layer: 1 (project-wide, no modality-specific logic)
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(yaml)
})

#' Load .env.sh variables into R session
load_env <- function(env_file = here(".env.sh")) {
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

#' Load dataset-specific config from configs_private
#' @param dataset_id character (optional, auto-detect from .env.sh)
#' @return list (parsed YAML)
load_dataset_config <- function(dataset_id = NULL) {
  if (is.null(dataset_id)) {
    dataset_id <- Sys.getenv("DATASET_ID", unset = "")
    if (!nzchar(dataset_id)) {
      load_env()
      dataset_id <- Sys.getenv("DATASET_ID", unset = "")
    }
  }
  if (!nzchar(dataset_id)) stop("DATASET_ID not set")
  
  cfg_path <- here("configs_private", "datasets", paste0(dataset_id, ".yaml"))
  if (!file.exists(cfg_path)) stop("Config not found: ", cfg_path)
  
  log_msg(sprintf("Loaded dataset config: %s", dataset_id))
  yaml::read_yaml(cfg_path)
}

#' Get DS_RAW path
get_ds_raw <- function(dataset_id = NULL) {
  # 1. Environment variable
  ds_raw <- Sys.getenv("DS_RAW", unset = "")
  if (nzchar(ds_raw) && dir.exists(ds_raw)) return(ds_raw)
  
  # 2. Load from .env.sh
  load_env()
  ds_raw <- Sys.getenv("DS_RAW", unset = "")
  if (nzchar(ds_raw) && dir.exists(ds_raw)) return(ds_raw)
  
  # 3. Private config
  cfg <- load_dataset_config(dataset_id)
  ds_raw <- cfg$project_data_root %||% ""
  if (nzchar(ds_raw) && dir.exists(ds_raw)) return(ds_raw)
  
  stop("DS_RAW not found. Source .env.sh or check configs_private/")
}

#' Standardized timestamp for filenames
timestamp_str <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

#' Logging helper
log_msg <- function(msg, level = "info") {
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%H:%M:%S"), level, msg))
}

cat("[init] scripts/utils/utils_io.R loaded\n")
