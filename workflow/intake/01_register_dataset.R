# workflow/intake/01_register_dataset.R
# Why:
# Read one dataset YAML config and upsert one machine-readable record into dataset_registry.csv.

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(tibble)
})

source(here("scripts", "utils_registry.R"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript workflow/intake/01_register_dataset.R <config.yaml>")
}

config_path <- args[[1]]

if (!file.exists(config_path)) {
  stop("Config file does not exist: ", config_path)
}

cfg <- yaml::read_yaml(config_path)

required_fields <- c(
  "dataset_id", "display_name", "dataset_type", "disease", "subtype",
  "modality", "data_level", "source_name", "owner", "status",
  "intake_qc_status", "harmonization_status", "analysis_qc_status",
  "atlas_role", "inclusion_decision", "download_date", "register_date",
  "last_update", "paths"
)

missing_fields <- setdiff(required_fields, names(cfg))
if (length(missing_fields) > 0) {
  stop("Missing required config fields: ", paste(missing_fields, collapse = ", "))
}

entry <- list(
  dataset_id = cfg$dataset_id,
  display_name = cfg$display_name,
  dataset_type = cfg$dataset_type,
  disease = cfg$disease,
  subtype = cfg$subtype,
  modality = cfg$modality,
  data_level = cfg$data_level,
  source_name = cfg$source_name,
  n_cohorts_claimed = cfg$n_cohorts_claimed,
  n_samples_claimed = cfg$n_samples_claimed,
  n_samples_aligned = cfg$n_samples_aligned,
  hpc_path = cfg$paths$raw_root,
  status = cfg$status,
  intake_qc_status = cfg$intake_qc_status,
  harmonization_status = cfg$harmonization_status,
  analysis_qc_status = cfg$analysis_qc_status,
  atlas_role = cfg$atlas_role,
  inclusion_decision = cfg$inclusion_decision,
  owner = cfg$owner,
  download_date = cfg$download_date,
  register_date = cfg$register_date,
  last_update = cfg$last_update,
  notes = cfg$notes
)

upsert_registry_entry(entry)
message("[OK] Registry upserted for: ", cfg$dataset_id)