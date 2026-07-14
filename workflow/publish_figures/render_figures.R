#!/usr/bin/env Rscript

command <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command, value = TRUE)
script_path <- normalizePath(
  sub("^--file=", "", file_arg[1]),
  winslash = "/",
  mustWork = TRUE
)
module_root <- dirname(script_path)
source(file.path(module_root, "R", "contracts.R"))
source(file.path(module_root, "R", "exporters.R"))
source(file.path(module_root, "R", "palettes.R"))
source(file.path(module_root, "R", "theme_cns.R"))
source(file.path(module_root, "R", "builders_atlas.R"))
source(file.path(module_root, "R", "builders_dotplot.R"))
source(file.path(module_root, "R", "builders_composition.R"))
source(file.path(module_root, "R", "render.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 4L) {
  stop(
    paste(
      "Usage: render_figures.R <figure_spec.yaml> <source_contracts.yaml>",
      "[project_root] [output_root]"
    ),
    call. = FALSE
  )
}

figure_spec <- load_figure_spec(normalizePath(args[[1]], mustWork = TRUE))
contract_registry <- load_contract_registry(normalizePath(args[[2]], mustWork = TRUE))
project_root <- if (length(args) >= 3L) args[[3]] else "."
output_root <- if (length(args) >= 4L) args[[4]] else NULL
manifest <- render_figure_spec(figure_spec, contract_registry, project_root, output_root)
cat(sprintf(
  "Figure rendering complete: %d panels across %d figure packages.\n",
  nrow(manifest), length(unique(manifest$figure_id))
))
