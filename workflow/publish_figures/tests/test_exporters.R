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
source(file.path(module_root, "R", "exporters.R"))

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(yaml)
})

temporary_root <- tempfile("publish_figure_exporters_")
dir.create(file.path(temporary_root, "objects"), recursive = TRUE)

counts <- Matrix(
  c(
    5, 4, 0, 0, 1, 0,
    0, 1, 6, 5, 0, 1,
    2, 2, 1, 1, 3, 3,
    0, 0, 0, 4, 5, 6
  ),
  nrow = 4,
  byrow = TRUE,
  sparse = TRUE,
  dimnames = list(
    c("GeneA", "GeneB", "GeneC", "GeneD"),
    paste0("cell", seq_len(6))
  )
)
object <- CreateSeuratObject(counts)
object$sample_id <- c("S1", "S1", "S2", "S2", "S3", "S3")
object$group <- c("G0", "G0", "G1", "G1", "G1", "G1")
object$celltype <- c("TypeA", "TypeA", "TypeB", "TypeB", "TypeA", "TypeB")
object <- NormalizeData(object, verbose = FALSE)
embedding <- matrix(
  seq_len(12) / 10,
  ncol = 2,
  dimnames = list(colnames(object), c("UMAP_1", "UMAP_2"))
)
object[["umap"]] <- CreateDimReducObject(
  embeddings = embedding,
  key = "UMAP_",
  assay = "RNA"
)
saveRDS(object, file.path(temporary_root, "objects", "synthetic.rds"))

write_yaml(
  list(marker_panel = list(Identity = c("GeneA", "GeneB"), Shared = "GeneC")),
  file.path(temporary_root, "features.yaml")
)
spec <- list(
  schema_version = 1,
  dataset_id = "synthetic_dataset",
  output_root = "figure_data",
  objects = list(
    atlas = list(
      path = "objects/synthetic.rds",
      sample_column = "sample_id",
      group_column = "group",
      assay = "RNA",
      layer = "data",
      reduction = "umap"
    )
  ),
  views = list(
    list(
      id = "atlas_l2",
      object = "atlas",
      annotation_column = "celltype",
      annotation_level = "L2",
      annotation_order = c("TypeA", "TypeB"),
      sample_order = c("S1", "S2", "S3"),
      outputs = c("umap", "dotplot", "composition"),
      feature_source = list(path = "features.yaml", set = "marker_panel"),
      composition_denominators = c("view", "object")
    )
  )
)
validate_export_spec(spec)
manifest <- export_figure_data(spec, temporary_root)

view_root <- file.path(temporary_root, "figure_data", "atlas_l2")
umap <- read.csv(file.path(view_root, "umap_cells.csv"), check.names = FALSE)
dotplot <- read.csv(file.path(view_root, "dotplot_expression.csv"), check.names = FALSE)
composition <- read.csv(file.path(view_root, "sample_composition.csv"), check.names = FALSE)

stopifnot(nrow(manifest) == 3L)
stopifnot(nrow(umap) == 6L)
stopifnot(all(c(
  "cell_id", "sample_id", "group_key", "annotation_level",
  "annotation", "dim1", "dim2"
) %in% names(umap)))
stopifnot(nrow(dotplot) == 6L)
stopifnot(identical(unique(dotplot$feature), c("GeneA", "GeneB", "GeneC")))
stopifnot(all(dotplot$pct_expr >= 0 & dotplot$pct_expr <= 100))
stopifnot(nrow(composition) == 12L)
stopifnot(identical(sort(unique(composition$denominator_scope)), c("object", "view")))
stopifnot(any(composition$n_cells == 0L))
stopifnot(all(composition$group_key %in% c("G0", "G1")))
stopifnot(!any(composition$group_key %in% c("Control", "Treatment")))

bad_spec <- spec
bad_spec$views[[1]]$annotation_order <- "TypeA"
bad_error <- tryCatch(
  {
    export_figure_data(bad_spec, temporary_root, "bad_output")
    FALSE
  },
  error = function(e) grepl("annotation_order omits", conditionMessage(e))
)
stopifnot(bad_error)

unlink(temporary_root, recursive = TRUE)
cat("publish_figures exporter tests: PASS\n")
