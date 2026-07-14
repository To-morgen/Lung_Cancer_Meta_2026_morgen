#!/usr/bin/env Rscript

command <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command, value = TRUE)
script_path <- normalizePath(
  sub("^--file=", "", file_arg[1]),
  winslash = "/",
  mustWork = TRUE
)
test_dir <- dirname(script_path)
module_root <- normalizePath(file.path(test_dir, ".."), winslash = "/")
project_root <- normalizePath(file.path(module_root, "..", ".."), winslash = "/")

source(file.path(module_root, "R", "contracts.R"))
source(file.path(module_root, "R", "palettes.R"))
source(file.path(module_root, "R", "theme_cns.R"))

registry_path <- file.path(
  module_root, "examples", "source_contracts.synthetic.yaml"
)
registry <- load_contract_registry(registry_path)
audit <- audit_contract_registry(registry, project_root)

stopifnot(nrow(audit) == 3L)
stopifnot(all(audit$ready[audit$state == "present"]))
stopifnot(!any(audit$blocking))
stopifnot(audit$state[audit$id == "validation_umap_cells"] == "derivable")

bad_path_error <- tryCatch(
  {
    resolve_contract_path(project_root, "/tmp/not-portable.csv")
    FALSE
  },
  error = function(e) TRUE
)
stopifnot(bad_path_error)

missing_column_contract <- list(
  id = "missing_column",
  state = "present",
  required_for = "main",
  privacy = "public",
  source = list(
    type = "table",
    path = "workflow/publish_figures/examples/fixtures/qc_sample_metrics.csv",
    format = "csv",
    required_columns = c("sample_id", "column_that_does_not_exist")
  )
)
missing_column_audit <- audit_source_contract(
  missing_column_contract,
  project_root
)
stopifnot(!missing_column_audit$ready)
stopifnot(missing_column_audit$blocking)
stopifnot(grepl("column_that_does_not_exist", missing_column_audit$reason))

palette <- name_palette(c("Type_A", "Type_B"))
stopifnot(identical(names(palette), c("Type_A", "Type_B")))
stopifnot(figure_width("single") < figure_width("double"))

cat("publish_figures contract tests: PASS\n")
