# workflow/intake/02_manifest_check.R
# Why:
# Build machine-readable manifest and md5 inventory for one dataset based on YAML config.

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
})

source(here("scripts", "utils_registry.R"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript workflow/intake/02_manifest_check.R <config.yaml>")
}

config_path <- args[[1]]

if (!file.exists(config_path)) {
  stop("Config file does not exist: ", config_path)
}

cfg <- yaml::read_yaml(config_path)

raw_root <- here(cfg$paths$raw_root)
manifest_root <- here(cfg$paths$manifest_root)

dir.create(manifest_root, recursive = TRUE, showWarnings = FALSE)

file_keys <- names(cfg$files)
file_names <- unlist(cfg$files, use.names = FALSE)
file_paths <- file.path(raw_root, file_names)

manifest_df <- tibble(
  file_key = file_keys,
  file_name = file_names,
  file_path = file_paths,
  exists = file.exists(file_paths),
  size_bytes = ifelse(file.exists(file_paths), file.info(file_paths)$size, NA_real_),
  mtime = ifelse(file.exists(file_paths), as.character(file.info(file_paths)$mtime), NA_character_),
  md5 = ifelse(file.exists(file_paths), unname(tools::md5sum(file_paths)), NA_character_)
)

readr::write_csv(manifest_df, file.path(manifest_root, "file_manifest.csv"))

if (all(manifest_df$exists)) {
  writeLines(
    paste(manifest_df$md5, manifest_df$file_name),
    con = file.path(manifest_root, "md5.txt")
  )
} else {
  missing_files <- manifest_df$file_name[!manifest_df$exists]
  warning("Missing files detected: ", paste(missing_files, collapse = ", "))
}

update_registry_fields(
  target_dataset_id = cfg$dataset_id,
  fields = list(
    last_update = as.character(Sys.Date()),
    notes = if (all(manifest_df$exists)) {
      "Manifest check completed"
    } else {
      paste0("Manifest check completed; missing files: ", paste(manifest_df$file_name[!manifest_df$exists], collapse = ", "))
    }
  )
)

message("[OK] Manifest written to: ", manifest_root)