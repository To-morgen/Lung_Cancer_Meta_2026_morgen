#!/usr/bin/env Rscript

command <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
module_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
source(file.path(module_root, "R", "contracts.R"))
source(file.path(module_root, "R", "exporters.R"))
source(file.path(module_root, "R", "palettes.R"))
source(file.path(module_root, "R", "theme_cns.R"))
source(file.path(module_root, "R", "builders_atlas.R"))
source(file.path(module_root, "R", "builders_dotplot.R"))
source(file.path(module_root, "R", "builders_composition.R"))
source(file.path(module_root, "R", "render.R"))

temporary_root <- tempfile("publish_figure_builders_")
dir.create(file.path(temporary_root, "inputs"), recursive = TRUE)

umap <- data.frame(
  cell_id = paste0("cell", 1:8),
  sample_id = rep(c("S1", "S2", "S3", "S4"), each = 2),
  group_key = rep(c("G0", "G0", "G1", "G1"), each = 2),
  annotation_level = "L2",
  annotation = rep(c("TypeA", "TypeB"), 4),
  dim1 = c(-2, -1, -1.5, -0.5, 0.5, 1.5, 1, 2),
  dim2 = c(0, 1, -1, 0, 0, 1, -1, 0),
  stringsAsFactors = FALSE
)
dotplot <- expand.grid(
  annotation = c("TypeA", "TypeB"),
  feature = c("GeneA", "GeneB"),
  stringsAsFactors = FALSE
)
dotplot$annotation_level <- "L2"
dotplot$annotation_order <- match(dotplot$annotation, c("TypeA", "TypeB"))
dotplot$feature_group <- "Identity"
dotplot$feature_order <- match(dotplot$feature, c("GeneA", "GeneB"))
dotplot$avg_expr <- c(2, 1, 1, 2)
dotplot$avg_expr_scaled <- c(1, -1, -1, 1)
dotplot$pct_expr <- c(80, 20, 25, 75)
dotplot$n_cells <- 4
dotplot$assay <- "RNA"
dotplot$layer <- "data"
composition <- expand.grid(
  sample_id = c("S1", "S2", "S3", "S4"),
  annotation = c("TypeA", "TypeB"),
  stringsAsFactors = FALSE
)
composition$group_key <- ifelse(composition$sample_id %in% c("S1", "S2"), "G0", "G1")
composition$annotation_level <- "L2"
composition$n_cells <- c(60, 55, 30, 25, 40, 45, 70, 75)
composition$denominator_cells <- 100
composition$denominator_scope <- "object"
composition$proportion <- composition$n_cells / composition$denominator_cells

write.csv(umap, file.path(temporary_root, "inputs", "umap.csv"), row.names = FALSE)
write.csv(dotplot, file.path(temporary_root, "inputs", "dotplot.csv"), row.names = FALSE)
write.csv(composition, file.path(temporary_root, "inputs", "composition.csv"), row.names = FALSE)

registry <- list(
  schema_version = 1,
  contracts = list(
    list(
      id = "umap", state = "present", required_for = "main", privacy = "public",
      source = list(type = "table", path = "inputs/umap.csv", format = "csv")
    ),
    list(
      id = "dotplot", state = "present", required_for = "main", privacy = "public",
      source = list(type = "table", path = "inputs/dotplot.csv", format = "csv")
    ),
    list(
      id = "composition", state = "present", required_for = "main", privacy = "public",
      source = list(type = "table", path = "inputs/composition.csv", format = "csv")
    )
  )
)
figure_spec <- list(
  schema_version = 1,
  manuscript_id = "synthetic",
  dataset_id = "synthetic",
  output_root = "rendered",
  display_groups = list(G0 = "Reference", G1 = "Intervention"),
  figures = list(
    list(
      id = "DEMO",
      panels = list(
        list(
          id = "A", builder = "atlas_umap", source_contracts = "umap",
          parameters = list(
            annotation_order = c("TypeA", "TypeB"),
            annotation_labels = list(TypeA = "Type A", TypeB = "Type B"),
            annotation_palette = list(TypeA = "#0072B2", TypeB = "#D55E00"),
            label = TRUE
          ),
          dimensions = list(width_in = 3, height_in = 2.6),
          output_formats = "png"
        ),
        list(
          id = "B", builder = "annotation_dotplot", source_contracts = "dotplot",
          dimensions = list(width_in = 3.5, height_in = 2.6),
          output_formats = "png"
        ),
        list(
          id = "C", builder = "sample_composition", source_contracts = "composition",
          parameters = list(denominator_scope = "object", ncol = 2),
          dimensions = list(width_in = 4, height_in = 2.8),
          output_formats = "png"
        )
      )
    )
  )
)

manifest <- render_figure_spec(figure_spec, registry, temporary_root)
stopifnot(nrow(manifest) == 3L)
stopifnot(all(file.exists(file.path(
  temporary_root,
  "rendered",
  "DEMO",
  c("A_atlas_umap.png", "B_annotation_dotplot.png", "C_sample_composition.png")
))))
stopifnot(file.exists(file.path(temporary_root, "rendered", "render_manifest.csv")))

unlink(temporary_root, recursive = TRUE)
cat("publish_figures builder tests: PASS\n")
