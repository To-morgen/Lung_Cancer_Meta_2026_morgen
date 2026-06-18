#!/usr/bin/env Rscript
# ============================================================================
# 06_prep_infercnv_inputs.R — Prepare inferCNV inputs
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(yaml)
  library(data.table)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

as_character_vec <- function(x) {
  if (is.null(x)) character() else as.character(unlist(x, use.names = FALSE))
}

as_single_integer <- function(x, default = NA_integer_) {
  if (is.null(x) || length(x) == 0) return(default)
  value <- as.integer(x[[1]])
  if (is.na(value)) default else value
}

proj_root <- Sys.getenv("LUNGMETA_ROOT", unset = "")
if (!nzchar(proj_root)) {
  cwd <- normalizePath(getwd())
  if (dir.exists(file.path(cwd, "modules", "cnv"))) {
    proj_root <- cwd
  } else {
    proj_root <- normalizePath(file.path(cwd, "..", ".."))
  }
}
module_root <- file.path(proj_root, "modules", "cnv")

path_from_root <- function(path) {
  if (is.null(path) || !nzchar(path)) return(path)
  if (grepl("^/", path)) path else file.path(proj_root, path)
}

resolve_path <- function(path) {
  if (is.null(path) || !nzchar(path)) return(path)
  if (grepl("^/", path)) return(path)
  project_path <- file.path(proj_root, path)
  module_path <- file.path(module_root, path)
  if (file.exists(project_path) || dir.exists(dirname(project_path))) return(project_path)
  if (file.exists(module_path) || dir.exists(dirname(module_path))) return(module_path)
  project_path
}

cap_cells_by_group <- function(cells, group_values, max_per_group, seed, label) {
  max_per_group <- as_single_integer(max_per_group)
  if (is.na(max_per_group)) return(cells)
  if (max_per_group < 1) stop("Invalid ", label, " cap: ", max_per_group)

  cells <- unique(cells)
  groups <- group_values[cells]
  if (any(is.na(groups))) stop("Missing group values while capping ", label)

  set.seed(seed)
  kept <- unlist(lapply(split(cells, groups), function(group_cells) {
    group_cells <- sort(group_cells)
    if (length(group_cells) > max_per_group) sample(group_cells, max_per_group) else group_cells
  }), use.names = FALSE)

  sort(kept)
}

get_counts <- function(sobj) {
  tryCatch(
    GetAssayData(sobj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(sobj, assay = "RNA", slot = "counts")
  )
}

named_counts_list <- function(values) {
  tbl <- table(values)
  out <- as.list(as.integer(tbl))
  names(out) <- names(tbl)
  out
}

cat("========================================\n")
cat("inferCNV Input Preparation\n")
cat(sprintf("Time: %s\n", Sys.time()))
cat("========================================\n")

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else if (nzchar(Sys.getenv("CNV_CONFIG", unset = ""))) {
  Sys.getenv("CNV_CONFIG")
} else {
  file.path(module_root, "configs", "cnv_params.yaml")
}
if (!file.exists(cfg_path)) stop("Config not found: ", cfg_path)

cfg <- yaml::read_yaml(cfg_path)
dataset_id <- cfg$dataset_id %||% "unspecified"
species <- cfg$species %||% "unspecified"
seed <- as_single_integer(cfg$infercnv$random_seed, 42L)

cat(sprintf("Config: %s\n", cfg_path))
cat(sprintf("Dataset: %s | Species: %s\n", dataset_id, species))

is_cnv_target_config <- !is.null(cfg$input$seurat_object) && !is.null(cfg$infercnv$observation_clusters)
is_legacy_role_config <- !is.null(cfg$seurat_object) && !is.null(cfg$role_column)

if (!is_cnv_target_config && !is_legacy_role_config) {
  stop("Config must be a CNV target YAML or legacy cnv_role dataset YAML")
}

if (is_cnv_target_config) {
  input_dir <- path_from_root(cfg$output$infercnv_inputs %||% file.path(cfg$output$base, "infercnv", "inputs"))
  seurat_path <- path_from_root(cfg$input$seurat_object)
  cluster_col <- cfg$input$cluster_column %||% "seurat_clusters"
  l2_col <- cfg$input$celltype_l2_column
  if (is.null(l2_col)) stop("Missing input.celltype_l2_column in config")

  ref_clusters <- as_character_vec(cfg$infercnv$reference_clusters)
  obs_clusters <- as_character_vec(cfg$infercnv$observation_clusters)
  ref_l2_include <- as_character_vec(cfg$infercnv$reference_l2_include)
  if (length(ref_clusters) == 0) stop("infercnv.reference_clusters is empty")
  if (length(obs_clusters) == 0) stop("infercnv.observation_clusters is empty")
  if (length(ref_l2_include) == 0) stop("infercnv.reference_l2_include is empty")

  overlap_clusters <- intersect(ref_clusters, obs_clusters)
  if (length(overlap_clusters) > 0) {
    stop("Reference and observation clusters overlap: ", paste(overlap_clusters, collapse = ","))
  }

  gene_order_path <- resolve_path(cfg$infercnv$gene_order_file %||% cfg$gene_order_file)
  mode <- "cnv_target_yaml"
} else {
  input_dir <- resolve_path(cfg$output$infercnv_inputs)
  seurat_path <- resolve_path(cfg$seurat_object)
  role_col <- cfg$role_column
  group_col <- cfg$group_column
  gene_order_path <- resolve_path(cfg$gene_order_file)
  mode <- "legacy_cnv_role"
}

if (is.null(input_dir) || !nzchar(input_dir)) stop("Missing inferCNV input output directory")
if (is.null(seurat_path) || !file.exists(seurat_path)) stop("Seurat object not found: ", seurat_path)
if (is.null(gene_order_path) || !file.exists(gene_order_path)) stop("Gene order not found: ", gene_order_path)

dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Mode: %s\n", mode))
cat(sprintf("Loading: %s\n", seurat_path))
sobj <- readRDS(seurat_path)
tryCatch(
  sobj <- JoinLayers(sobj),
  error = function(e) cat(sprintf("  JoinLayers skipped: %s\n", e$message))
)
cat(sprintf("  Cells: %d | Genes: %d\n", ncol(sobj), nrow(sobj)))

md <- sobj@meta.data

if (is_cnv_target_config) {
  if (!cluster_col %in% colnames(md)) stop("Missing metadata column: ", cluster_col)
  if (!l2_col %in% colnames(md)) stop("Missing metadata column: ", l2_col)

  clusters <- setNames(as.character(md[[cluster_col]]), rownames(md))
  l2_values <- setNames(as.character(md[[l2_col]]), rownames(md))

  cells_ref <- rownames(md)[clusters %in% ref_clusters & l2_values %in% ref_l2_include]
  cells_obs <- rownames(md)[clusters %in% obs_clusters]

  if (length(cells_ref) == 0) stop("No reference cells selected")
  if (length(cells_obs) == 0) stop("No observation cells selected")

  overlap_cells <- intersect(cells_ref, cells_obs)
  if (length(overlap_cells) > 0) stop("Reference and observation cells overlap: ", length(overlap_cells))

  cells_ref <- cap_cells_by_group(
    cells_ref,
    l2_values,
    cfg$infercnv$max_cells_per_reference_l2,
    seed,
    "reference L2"
  )
  cells_obs <- cap_cells_by_group(
    cells_obs,
    clusters,
    cfg$infercnv$max_cells_per_observation_cluster,
    seed,
    "observation cluster"
  )

  group_map <- character(length(cells_ref) + length(cells_obs))
  names(group_map) <- c(cells_ref, cells_obs)
  group_map[cells_ref] <- paste0("ref_", l2_values[cells_ref])
  group_map[cells_obs] <- paste0("obs_C", clusters[cells_obs])
  cells_excl <- setdiff(rownames(md), names(group_map))
} else {
  if (!role_col %in% colnames(md)) stop("Missing metadata column: ", role_col)
  if (!group_col %in% colnames(md)) stop("Missing metadata column: ", group_col)
  if (role_col != "cnv_role") stop("Legacy mode expects role_column: cnv_role")

  cells_ref <- rownames(md[md[[role_col]] == "reference", ])
  cells_obs <- rownames(md[md[[role_col]] == "observation", ])
  cells_excl <- rownames(md[grepl("excluded", md[[role_col]]), ])

  if (length(cells_ref) == 0) stop("No reference cells selected")
  if (length(cells_obs) == 0) stop("No observation cells selected")

  group_values <- setNames(as.character(md[[group_col]]), rownames(md))
  group_map <- group_values[c(cells_ref, cells_obs)]
}

cat(sprintf("  Reference:   %d\n  Observation: %d\n  Excluded:    %d\n", length(cells_ref), length(cells_obs), length(cells_excl)))

cells_keep <- unique(c(cells_ref, cells_obs))
sobj <- subset(sobj, cells = cells_keep)
md <- sobj@meta.data

counts <- get_counts(sobj)
cat(sprintf("  Counts: %d x %d\n", nrow(counts), ncol(counts)))

annotations <- data.frame(
  barcode = colnames(sobj),
  group = unname(group_map[colnames(sobj)]),
  stringsAsFactors = FALSE
)
if (any(is.na(annotations$group))) stop("Missing inferCNV group assignments after subsetting")

ref_groups <- sort(unique(annotations$group[annotations$barcode %in% cells_ref]))
obs_groups <- sort(unique(annotations$group[annotations$barcode %in% cells_obs]))
if (length(ref_groups) == 0) stop("No reference groups in annotation")
if (length(obs_groups) == 0) stop("No observation groups in annotation")

cat(sprintf("\n  Ref groups (%d):\n", length(ref_groups)))
for (g in ref_groups) cat(sprintf("    %s: %d\n", g, sum(annotations$group == g)))
cat(sprintf("  Obs groups (%d):\n", length(obs_groups)))
for (g in obs_groups) cat(sprintf("    %s: %d\n", g, sum(annotations$group == g)))

gene_order <- fread(gene_order_path, header = FALSE)
dup_count <- sum(duplicated(gene_order$V1))
cat(sprintf("\n  Gene order: %d | Duplicates: %d\n", nrow(gene_order), dup_count))
if (dup_count > 0) stop("Gene order file contains duplicated gene symbols: ", dup_count)

overlap <- intersect(rownames(counts), gene_order$V1)
cat(sprintf("  Overlap: %d (%.1f%%)\n", length(overlap), 100 * length(overlap) / nrow(counts)))
if (length(overlap) == 0) stop("No overlap between counts genes and gene order")

saveRDS(counts, file.path(input_dir, "counts_sparse.rds"))
write.table(
  annotations,
  file.path(input_dir, "cell_annotations.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
writeLines(ref_groups, file.path(input_dir, "reference_groups.txt"))

yaml::write_yaml(list(
  dataset_id = dataset_id,
  species = species,
  mode = mode,
  source_config = normalizePath(cfg_path),
  seurat_object = normalizePath(seurat_path),
  gene_order_file = normalizePath(gene_order_path),
  n_ref_cells = length(cells_ref),
  n_obs_cells = length(cells_obs),
  n_excluded = length(cells_excl),
  n_genes = nrow(counts),
  gene_overlap = length(overlap),
  ref_groups = ref_groups,
  obs_groups = obs_groups,
  ref_group_counts = named_counts_list(annotations$group[annotations$barcode %in% cells_ref]),
  obs_group_counts = named_counts_list(annotations$group[annotations$barcode %in% cells_obs])
), file.path(input_dir, "dataset_info.yaml"))

fwrite(data.frame(
  metric = c("dataset_id", "species", "mode", "total_cells", "total_genes", "n_ref", "n_obs", "n_excl", "gene_overlap"),
  value = c(dataset_id, species, mode, ncol(counts), nrow(counts), length(cells_ref), length(cells_obs), length(cells_excl), length(overlap))
), file.path(input_dir, "prep_summary.csv"))

cat("\nDone: ", dataset_id, " -> ", input_dir, "\n", sep = "")
