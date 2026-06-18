#!/usr/bin/env Rscript
# ============================================================================
# 02_infercnv.R — Run inferCNV from prepared inputs
# ============================================================================

suppressPackageStartupMessages({
  library(infercnv)
  library(Matrix)
  library(data.table)
  library(yaml)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

as_single_integer <- function(x, default = NA_integer_) {
  if (is.null(x) || length(x) == 0) return(default)
  value <- as.integer(x[[1]])
  if (is.na(value)) default else value
}

as_single_logical <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0) return(default)
  value <- as.logical(x[[1]])
  if (is.na(value)) default else value
}

merge_lists <- function(defaults, overrides) {
  if (is.null(defaults)) defaults <- list()
  if (is.null(overrides)) return(defaults)
  for (name in names(overrides)) defaults[[name]] <- overrides[[name]]
  defaults
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

cat("\n")
cat("╔══════════════════════════════════════════════════╗\n")
cat("║   CNV Module: inferCNV Inference                 ║\n")
cat(sprintf("║   Time: %s            ║\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("╚══════════════════════════════════════════════════╝\n\n")

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
shared_path <- file.path(module_root, "configs", "cnv_params_shared.yaml")
shared_cfg <- if (file.exists(shared_path)) yaml::read_yaml(shared_path) else list(infercnv = list())

run_cfg <- merge_lists(shared_cfg$infercnv, cfg$infercnv)
selection_keys <- c(
  "enabled", "strategy", "reference_strategy", "reference_clusters",
  "reference_l2_include", "observation_clusters", "max_cells_per_reference_l2",
  "max_cells_per_observation_cluster", "gene_order_file", "random_seed", "container"
)
run_cfg[intersect(names(run_cfg), selection_keys)] <- NULL

threads_env <- Sys.getenv("INFERCNV_THREADS", unset = "")
if (nzchar(threads_env)) {
  run_cfg$num_threads <- as_single_integer(threads_env)
}
if (is.null(run_cfg$num_threads) || is.na(as.integer(run_cfg$num_threads)) || as.integer(run_cfg$num_threads) < 1) {
  run_cfg$num_threads <- 1L
}

input_dir <- path_from_root(cfg$output$infercnv_inputs %||% file.path(cfg$output$base, "infercnv", "inputs"))
output_dir <- path_from_root(cfg$output$infercnv_run %||% file.path(cfg$output$base, "infercnv", "run"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- file.path(input_dir, c(
  "counts_sparse.rds",
  "cell_annotations.tsv",
  "reference_groups.txt",
  "dataset_info.yaml"
))
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("Missing prepared inferCNV input(s):\n", paste(missing_inputs, collapse = "\n"))
}

ds_info <- yaml::read_yaml(file.path(input_dir, "dataset_info.yaml"))
gene_order <- ds_info$gene_order_file
if (is.null(gene_order) || !file.exists(gene_order)) stop("Gene order not found from dataset_info.yaml: ", gene_order)

cat(sprintf("Config:     %s\n", cfg_path))
cat(sprintf("Dataset:    %s\n", dataset_id))
cat(sprintf("Input dir:  %s\n", input_dir))
cat(sprintf("Output dir: %s\n", output_dir))
cat(sprintf("Gene order: %s\n", gene_order))
cat(sprintf("Threads:    %d\n", as.integer(run_cfg$num_threads)))
cat(sprintf("Mode:       %s\n", run_cfg$analysis_mode %||% "default"))
cat(sprintf("HMM:        %s\n", as.character(run_cfg$HMM %||% FALSE)))

cat("\nLoading prepared inputs...\n")
counts <- readRDS(file.path(input_dir, "counts_sparse.rds"))
ref_groups <- readLines(file.path(input_dir, "reference_groups.txt"))
annot_path <- file.path(input_dir, "cell_annotations.tsv")
cat(sprintf("  Counts: %d genes x %d cells\n", nrow(counts), ncol(counts)))
cat(sprintf("  Reference groups: %s\n", paste(ref_groups, collapse = ", ")))

cat("\nCreating inferCNV object...\n")
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file = annot_path,
  gene_order_file = gene_order,
  ref_group_names = ref_groups
)
cat(sprintf("  inferCNV object: %d genes x %d cells\n", nrow(infercnv_obj@expr.data), ncol(infercnv_obj@expr.data)))

cat("\n══════════════════════════════════════════════════\n")
cat(sprintf("   Running inferCNV [%s]\n", format(Sys.time(), "%H:%M:%S")))
cat("══════════════════════════════════════════════════\n\n")

t0 <- Sys.time()

infercnv_obj <- infercnv::run(
  infercnv_obj,
  min_cells_per_gene = run_cfg$min_cells_per_gene %||% NULL,
  cutoff = run_cfg$cutoff %||% 0.1,
  out_dir = output_dir,
  cluster_by_groups = as_single_logical(run_cfg$cluster_by_groups, TRUE),
  analysis_mode = run_cfg$analysis_mode %||% "subclusters",
  tumor_subcluster_partition_method = run_cfg$tumor_subcluster_partition_method %||% "leiden",
  tumor_subcluster_pval = run_cfg$tumor_subcluster_pval %||% 0.1,
  k_nn = run_cfg$k_nn %||% 20,
  denoise = as_single_logical(run_cfg$denoise, TRUE),
  noise_filter = run_cfg$noise_filter %||% 0.05,
  sd_amplifier = run_cfg$sd_amplifier %||% 1.0,
  HMM = as_single_logical(run_cfg$HMM, FALSE),
  HMM_type = run_cfg$HMM_type %||% "i6",
  BayesMaxPNormal = run_cfg$BayesMaxPNormal %||% 0.5,
  reassignCNVs = as_single_logical(run_cfg$reassignCNVs, TRUE),
  num_threads = as.integer(run_cfg$num_threads),
  output_format = run_cfg$output_format %||% "png",
  plot_steps = as_single_logical(run_cfg$plot_steps, FALSE),
  no_prelim_plot = as_single_logical(run_cfg$no_prelim_plot, FALSE),
  resume_mode = as_single_logical(run_cfg$resume_mode, TRUE),
  png_res = run_cfg$png_res %||% 300,
  no_plot = as_single_logical(run_cfg$no_plot, FALSE)
)

elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
cat(sprintf("\ninferCNV completed in %.1f min\n", elapsed))

saveRDS(infercnv_obj, file.path(output_dir, "infercnv_obj_final.rds"))
yaml::write_yaml(list(
  dataset_id = dataset_id,
  source_config = normalizePath(cfg_path),
  input_dir = normalizePath(input_dir),
  output_dir = normalizePath(output_dir),
  gene_order_file = normalizePath(gene_order),
  reference_groups = ref_groups,
  n_genes = nrow(infercnv_obj@expr.data),
  n_cells = ncol(infercnv_obj@expr.data),
  elapsed_minutes = as.numeric(elapsed),
  num_threads = as.integer(run_cfg$num_threads),
  analysis_mode = run_cfg$analysis_mode %||% "subclusters",
  HMM = as_single_logical(run_cfg$HMM, FALSE),
  HMM_type = run_cfg$HMM_type %||% "i6"
), file.path(output_dir, "run_summary.yaml"))

cat("\n")
cat("╔══════════════════════════════════════════════════════╗\n")
cat("║                inferCNV Complete                     ║\n")
cat(sprintf("║  Elapsed: %.1f min\n", elapsed))
cat(sprintf("║  Output:  %s\n", output_dir))
cat("╚══════════════════════════════════════════════════════╝\n")

gc()
