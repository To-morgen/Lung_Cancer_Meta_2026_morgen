#!/usr/bin/env Rscript

command <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command, value = TRUE)
script_path <- normalizePath(
  sub("^--file=", "", file_arg[1]),
  winslash = "/",
  mustWork = TRUE
)
module_root <- dirname(script_path)
source(file.path(module_root, "R", "exporters.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 3L) {
  stop(
    paste(
      "Usage: export_figure_data.R <export_spec.yaml>",
      "[project_root] [output_root]"
    ),
    call. = FALSE
  )
}

spec_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
project_root <- if (length(args) >= 2L) args[[2]] else "."
output_root <- if (length(args) >= 3L) args[[3]] else NULL

spec <- load_export_spec(spec_path)
manifest <- export_figure_data(spec, project_root, output_root)
cat(sprintf(
  "Figure-data export complete: %d outputs across %d views.\n",
  nrow(manifest), length(unique(manifest$view_id))
))
