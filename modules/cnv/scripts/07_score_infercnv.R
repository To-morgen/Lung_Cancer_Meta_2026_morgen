#!/usr/bin/env Rscript
# ============================================================================
# 07_score_infercnv.R — Score per-cell CNV burden from inferCNV output
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(yaml)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

as_character_vec <- function(x) {
  if (is.null(x)) character() else as.character(unlist(x, use.names = FALSE))
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

read_annotation <- function(path) {
  annot <- fread(path, header = FALSE, col.names = c("cell", "infercnv_group"))
  annot$cell <- as.character(annot$cell)
  annot$infercnv_group <- as.character(annot$infercnv_group)
  annot
}

infer_cluster <- function(group) {
  out <- sub("^obs_C", "", group)
  out[!grepl("^obs_C", group)] <- NA_character_
  out
}

score_from_expr <- function(mat) {
  mat_min <- suppressWarnings(min(mat, na.rm = TRUE))
  mat_max <- suppressWarnings(max(mat, na.rm = TRUE))
  if (is.finite(mat_min) && is.finite(mat_max) && mat_min >= 0 && mat_max <= 3) {
    colMeans((mat - 1)^2, na.rm = TRUE)
  } else {
    colMeans(abs(mat), na.rm = TRUE)
  }
}

summary_by <- function(df, group_cols, threshold_95, threshold_99) {
  df %>%
    filter(if_all(all_of(group_cols), ~ !is.na(.x))) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_cells = n(),
      mean_cnv_score = round(mean(cnv_score, na.rm = TRUE), 6),
      median_cnv_score = round(median(cnv_score, na.rm = TRUE), 6),
      q75_cnv_score = round(as.numeric(quantile(cnv_score, 0.75, na.rm = TRUE)), 6),
      q95_cnv_score = round(as.numeric(quantile(cnv_score, 0.95, na.rm = TRUE)), 6),
      n_above_ref95 = sum(cnv_score >= threshold_95, na.rm = TRUE),
      pct_above_ref95 = round(n_above_ref95 / n() * 100, 1),
      n_above_ref99 = sum(cnv_score >= threshold_99, na.rm = TRUE),
      pct_above_ref99 = round(n_above_ref99 / n() * 100, 1),
      .groups = "drop"
    ) %>%
    arrange(desc(median_cnv_score))
}

plot_violin <- function(df, x_col, title, threshold_95, threshold_99) {
  ggplot(df, aes(x = reorder(.data[[x_col]], cnv_score, median, na.rm = TRUE), y = cnv_score, fill = cell_class)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
    geom_boxplot(width = 0.12, outlier.size = 0.2, alpha = 0.7) +
    geom_hline(yintercept = threshold_95, linetype = "dashed", color = "red3") +
    geom_hline(yintercept = threshold_99, linetype = "dotted", color = "red4") +
    coord_flip() +
    scale_fill_manual(values = c("reference" = "#4DAF4A", "observation" = "#E41A1C")) +
    labs(title = title, x = "", y = "CNV score", fill = "") +
    theme_minimal() +
    theme(legend.position = "top")
}

cat("\n")
cat("==============================================================\n")
cat("   inferCNV CNV score postprocess                             \n")
cat("==============================================================\n\n")

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
base_out <- path_from_root(cfg$output$base %||% file.path("results", "scrna", dataset_id, "07_cnv"))
infercnv_run <- path_from_root(cfg$output$infercnv_run %||% file.path(base_out, "infercnv", "run"))
infercnv_inputs <- path_from_root(cfg$output$infercnv_inputs %||% file.path(base_out, "infercnv", "inputs"))
score_out <- path_from_root(cfg$output$infercnv_score %||% file.path(base_out, "infercnv", "scoring"))
report_dir <- file.path(score_out, "reports")
plot_dir <- file.path(score_out, "plots")
for (d in c(score_out, report_dir, plot_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

infercnv_rds <- file.path(infercnv_run, "infercnv_obj_final.rds")
annot_path <- file.path(infercnv_inputs, "cell_annotations.tsv")
seurat_path <- path_from_root(cfg$input$seurat_object)
if (!file.exists(infercnv_rds)) stop("inferCNV object not found: ", infercnv_rds)
if (!file.exists(annot_path)) stop("Prepared annotation not found: ", annot_path)
if (!file.exists(seurat_path)) stop("Seurat object not found: ", seurat_path)

cat(sprintf("Config:     %s\n", cfg_path))
cat(sprintf("Dataset:    %s\n", dataset_id))
cat(sprintf("inferCNV:   %s\n", infercnv_rds))
cat(sprintf("Output:     %s\n", score_out))

infercnv_obj <- readRDS(infercnv_rds)
expr_mat <- infercnv_obj@expr.data
cat(sprintf("Matrix:     %d genes x %d cells\n", nrow(expr_mat), ncol(expr_mat)))

scores <- score_from_expr(expr_mat)
score_df <- data.frame(
  cell = names(scores),
  cnv_score = as.numeric(scores),
  stringsAsFactors = FALSE
)

annot <- read_annotation(annot_path)
score_df <- score_df %>% left_join(annot, by = "cell")
score_df$cell_class <- ifelse(grepl("^ref_", score_df$infercnv_group), "reference", "observation")
score_df$infercnv_cluster <- infer_cluster(score_df$infercnv_group)

ref_scores <- score_df$cnv_score[score_df$cell_class == "reference"]
obs_scores <- score_df$cnv_score[score_df$cell_class == "observation"]
threshold_95 <- as.numeric(quantile(ref_scores, 0.95, na.rm = TRUE))
threshold_99 <- as.numeric(quantile(ref_scores, 0.99, na.rm = TRUE))
score_df$cnv_score_class <- case_when(
  score_df$cell_class == "reference" ~ "reference",
  score_df$cnv_score >= threshold_99 ~ "obs_above_ref99",
  score_df$cnv_score >= threshold_95 ~ "obs_above_ref95",
  TRUE ~ "obs_below_ref95"
)

cat(sprintf("Reference score: median=%.6f, 95th=%.6f, 99th=%.6f\n", median(ref_scores), threshold_95, threshold_99))
cat(sprintf("Observation score: median=%.6f, mean=%.6f\n", median(obs_scores), mean(obs_scores)))

cat("Loading Seurat metadata...\n")
sobj <- readRDS(seurat_path)
md <- sobj@meta.data
md$cell <- rownames(md)
metadata_cols <- intersect(
  c("cell", cfg$input$cluster_column %||% "seurat_clusters", cfg$input$celltype_l1_column %||% "celltype_L1", cfg$input$celltype_l2_column %||% "celltype_L2", "sample_id", "group"),
  colnames(md)
)
score_df <- score_df %>% left_join(md[, metadata_cols, drop = FALSE], by = "cell")

cluster_col <- cfg$input$cluster_column %||% "seurat_clusters"
l1_col <- cfg$input$celltype_l1_column %||% "celltype_L1"
l2_col <- cfg$input$celltype_l2_column %||% "celltype_L2"
if (!cluster_col %in% colnames(score_df)) score_df[[cluster_col]] <- score_df$infercnv_cluster
score_df$cluster_label <- paste0("C", score_df[[cluster_col]])
score_df$target_label <- ifelse(is.na(score_df[[l2_col]]), score_df$infercnv_group, score_df[[l2_col]])

fwrite(score_df, file.path(report_dir, "infercnv_per_cell_scores.csv"))
fwrite(data.frame(
  metric = c("ref_median", "ref_95th", "ref_99th", "obs_median", "obs_mean", "n_reference", "n_observation"),
  value = c(median(ref_scores), threshold_95, threshold_99, median(obs_scores), mean(obs_scores), length(ref_scores), length(obs_scores))
), file.path(report_dir, "infercnv_score_thresholds.csv"))

by_group <- summary_by(score_df, c("infercnv_group"), threshold_95, threshold_99)
by_cluster <- summary_by(score_df, c("cluster_label"), threshold_95, threshold_99)
fwrite(by_group, file.path(report_dir, "infercnv_score_by_group.csv"))
fwrite(by_cluster, file.path(report_dir, "infercnv_score_by_cluster.csv"))

if (all(c(l1_col, l2_col) %in% colnames(score_df))) {
  by_l2 <- summary_by(score_df, c(l1_col, l2_col), threshold_95, threshold_99)
  fwrite(by_l2, file.path(report_dir, "infercnv_score_by_l2.csv"))
}
if (all(c("sample_id", "group", "cluster_label") %in% colnames(score_df))) {
  by_sample_cluster <- summary_by(score_df, c("group", "sample_id", "cluster_label"), threshold_95, threshold_99)
  fwrite(by_sample_cluster, file.path(report_dir, "infercnv_score_by_sample_cluster.csv"))
}

focus_clusters <- paste0("C", as_character_vec(cfg$infercnv$observation_clusters))
focus_clusters <- intersect(focus_clusters, unique(score_df$cluster_label))
focus_df <- score_df %>% filter(cluster_label %in% focus_clusters | cell_class == "reference")

p_density <- ggplot(score_df, aes(x = cnv_score, fill = cell_class)) +
  geom_density(alpha = 0.45) +
  geom_vline(xintercept = threshold_95, linetype = "dashed", color = "red3") +
  geom_vline(xintercept = threshold_99, linetype = "dotted", color = "red4") +
  scale_fill_manual(values = c("reference" = "#4DAF4A", "observation" = "#E41A1C")) +
  labs(title = "inferCNV score distribution", x = "CNV score", y = "Density", fill = "") +
  theme_minimal() +
  theme(legend.position = "top")

p_cluster <- plot_violin(
  focus_df %>% filter(cell_class == "observation"),
  "cluster_label",
  "CNV score by observation cluster",
  threshold_95,
  threshold_99
)

if (l2_col %in% colnames(score_df)) {
  p_l2 <- plot_violin(
    focus_df %>% filter(cell_class == "observation"),
    l2_col,
    "CNV score by L2 label",
    threshold_95,
    threshold_99
  )
} else {
  p_l2 <- ggplot() + theme_void()
}

if (all(c("sample_id", "group") %in% colnames(score_df))) {
  p_sample <- focus_df %>%
    filter(cell_class == "observation") %>%
    ggplot(aes(x = cluster_label, y = cnv_score, fill = group)) +
    geom_boxplot(outlier.size = 0.2) +
    geom_hline(yintercept = threshold_95, linetype = "dashed", color = "red3") +
    facet_wrap(~sample_id, scales = "free_x") +
    labs(title = "CNV score by cluster and sample", x = "", y = "CNV score", fill = "Group") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
} else {
  p_sample <- ggplot() + theme_void()
}

pdf(file.path(plot_dir, "01_infercnv_score_overview.pdf"), width = 18, height = 18)
print((p_density | p_cluster) / (p_l2 | p_sample))
dev.off()

pdf(file.path(plot_dir, "02_infercnv_score_by_cluster.pdf"), width = 10, height = 7)
print(p_cluster)
dev.off()

pdf(file.path(plot_dir, "03_infercnv_score_by_l2.pdf"), width = 12, height = 8)
print(p_l2)
dev.off()

pdf(file.path(plot_dir, "04_infercnv_score_by_sample_cluster.pdf"), width = 14, height = 10)
print(p_sample)
dev.off()

saveRDS(score_df, file.path(score_out, "infercnv_score_df.rds"))
yaml::write_yaml(list(
  dataset_id = dataset_id,
  source_config = normalizePath(cfg_path),
  infercnv_object = normalizePath(infercnv_rds),
  n_cells = nrow(score_df),
  n_reference = length(ref_scores),
  n_observation = length(obs_scores),
  score_formula = "colMeans((expr - 1)^2) for nonnegative inferCNV final matrix; otherwise colMeans(abs(expr))",
  ref_median = as.numeric(median(ref_scores)),
  ref_95th = threshold_95,
  ref_99th = threshold_99,
  obs_median = as.numeric(median(obs_scores)),
  obs_mean = as.numeric(mean(obs_scores))
), file.path(score_out, "score_summary.yaml"))

cat("\nTop cluster summaries by median CNV score:\n")
print(as.data.frame(by_cluster %>% slice_head(n = 12)), row.names = FALSE)
cat("\nDone: ", score_out, "\n", sep = "")
