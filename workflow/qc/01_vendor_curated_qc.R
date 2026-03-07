# workflow/qc/01_vendor_curated_qc.R
# Why:
# Intake QC for vendor-curated bulk dataset:
# - read-check vendor objects
# - compare cohort coverage across cell / expr / clin layers
# - summarize object structure
# - stage bundle for downstream harmonization
# - update registry status

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(qs)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(stringr)
})

source(here("scripts", "utils_registry.R"))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript workflow/qc/01_vendor_curated_qc.R <config.yaml>")
}

config_path <- args[[1]]

if (!file.exists(config_path)) {
  stop("Config file does not exist: ", config_path)
}

cfg <- yaml::read_yaml(config_path)

raw_root <- here(cfg$paths$raw_root)
qc_root <- here(cfg$paths$qc_root)
stage_root <- here(cfg$paths$stage_root)

dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)
dir.create(stage_root, recursive = TRUE, showWarnings = FALSE)

safe_dim <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
    return(dim(x))
  }
  c(NA_integer_, NA_integer_)
}

guess_id_col <- function(df) {
  if (!is.data.frame(df)) {
    return(NA_character_)
  }

  cn <- colnames(df)
  low <- tolower(cn)

  patterns <- c("^sample_id$", "^sample$", "sample", "barcode", "patient", "^id$")
  for (pat in patterns) {
    hit <- cn[str_detect(low, pat)]
    if (length(hit) > 0) {
      return(hit[1])
    }
  }

  NA_character_
}

summarise_named_list <- function(x_list, object_type) {
  purrr::imap_dfr(x_list, function(x, nm) {
    obj_df <- if (is.matrix(x)) as.data.frame(x) else x
    d <- safe_dim(x)

    tibble(
      cohort = nm,
      object_type = object_type,
      class = class(x)[1],
      nrow = d[1],
      ncol = d[2],
      guessed_id_col = if (is.data.frame(obj_df)) guess_id_col(obj_df) else NA_character_
    )
  })
}

cell_qs_path <- file.path(raw_root, cfg$files$cell_qs)
symbol_qs_path <- file.path(raw_root, cfg$files$symbol_qs)

if (!file.exists(cell_qs_path)) {
  stop("Missing required file: ", cell_qs_path)
}

if (!file.exists(symbol_qs_path)) {
  stop("Missing required file: ", symbol_qs_path)
}

update_registry_fields(
  target_dataset_id = cfg$dataset_id,
  fields = list(
    status = "intake_qc",
    intake_qc_status = "pending",
    last_update = as.character(Sys.Date())
  )
)

cell_infilt_list <- qs::qread(cell_qs_path)
symbol_data <- qs::qread(symbol_qs_path)

if (!is.list(cell_infilt_list)) {
  stop("Cell infiltration object is not a list.")
}

if (!is.list(symbol_data) || !all(c("total_expr_list", "total_clin_list") %in% names(symbol_data))) {
  stop("symbol.qs does not contain expected keys: total_expr_list / total_clin_list")
}

expr_list <- symbol_data$total_expr_list
clin_list <- symbol_data$total_clin_list

coh_cell <- names(cell_infilt_list)
coh_expr <- names(expr_list)
coh_clin <- names(clin_list)

all_cohorts <- sort(unique(c(coh_cell, coh_expr, coh_clin)))

presence_map <- tibble(cohort = all_cohorts) %>%
  mutate(
    has_cell = cohort %in% coh_cell,
    has_expr = cohort %in% coh_expr,
    has_clin = cohort %in% coh_clin,
    presence_status = case_when(
      has_cell & has_expr & has_clin ~ "all_present",
      !has_cell & has_expr & has_clin ~ "missing_cell",
      has_cell & !has_expr & has_clin ~ "missing_expr",
      has_cell & has_expr & !has_clin ~ "missing_clin",
      TRUE ~ "other_mismatch"
    )
  )

object_summary <- bind_rows(
  summarise_named_list(cell_infilt_list, "cell"),
  summarise_named_list(expr_list, "expr"),
  summarise_named_list(clin_list, "clin")
) %>%
  arrange(cohort, object_type)

n_missing_cell <- sum(presence_map$presence_status == "missing_cell")
n_missing_expr <- sum(presence_map$presence_status == "missing_expr")
n_missing_clin <- sum(presence_map$presence_status == "missing_clin")

qc_status <- case_when(
  n_missing_expr > 0 | n_missing_clin > 0 ~ "fail",
  n_missing_cell > 0 ~ "warning",
  TRUE ~ "pass"
)

qc_note <- case_when(
  qc_status == "fail" ~ paste0(
    "Critical structural mismatch: missing_expr=", n_missing_expr,
    "; missing_clin=", n_missing_clin
  ),
  qc_status == "warning" ~ paste0(
    "Vendor object mismatch: missing_cell=", n_missing_cell,
    "; proceed conditionally to harmonization"
  ),
  TRUE ~ "All structural layers are consistent"
)

intake_summary <- tibble(
  dataset_id = cfg$dataset_id,
  n_cell_cohorts = length(coh_cell),
  n_expr_cohorts = length(coh_expr),
  n_clin_cohorts = length(coh_clin),
  n_all_unique_cohorts = length(all_cohorts),
  n_missing_cell = n_missing_cell,
  n_missing_expr = n_missing_expr,
  n_missing_clin = n_missing_clin,
  intake_qc_status = qc_status,
  note = qc_note
)

readr::write_csv(presence_map, file.path(qc_root, "cohort_presence_map.csv"))
readr::write_csv(object_summary, file.path(qc_root, "cohort_object_summary.csv"))
readr::write_csv(intake_summary, file.path(qc_root, "intake_qc_summary.csv"))

report_lines <- c(
  paste0("# Intake QC Report: ", cfg$dataset_id),
  "",
  paste0("- Date: ", Sys.Date()),
  paste0("- Dataset: ", cfg$display_name),
  paste0("- Cell cohorts: ", length(coh_cell)),
  paste0("- Expr cohorts: ", length(coh_expr)),
  paste0("- Clin cohorts: ", length(coh_clin)),
  paste0("- QC status: ", qc_status),
  paste0("- Note: ", qc_note)
)

writeLines(report_lines, con = file.path(qc_root, "intake_qc_report.md"))

saveRDS(
  list(
    dataset_id = cfg$dataset_id,
    config = cfg,
    cell_infilt_list = cell_infilt_list,
    expr_list = expr_list,
    clin_list = clin_list,
    presence_map = presence_map,
    object_summary = object_summary,
    intake_summary = intake_summary
  ),
  file = file.path(stage_root, "intake_readcheck_bundle.rds")
)

update_registry_fields(
  target_dataset_id = cfg$dataset_id,
  fields = list(
    status = "intake_qc",
    intake_qc_status = qc_status,
    last_update = as.character(Sys.Date()),
    notes = paste0("Vendor intake QC ", qc_status, ": ", qc_note)
  )
)

message("[OK] Vendor QC completed for: ", cfg$dataset_id)
message("[QC status] ", qc_status)