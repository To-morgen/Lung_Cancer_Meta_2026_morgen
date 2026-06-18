#!/usr/bin/env Rscript
# ============================================================================
# 05_summarize_reports.R -- Summarize generated CellChat CSV reports
# ============================================================================

cat(
  "\n+================================================================+\n",
  "|             CellChat Step 05: Summarize Reports              |\n",
  "+================================================================+\n",
  sep = ""
)

suppressPackageStartupMessages({
  library(yaml)
  library(here)
  library(data.table)
})

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

group_label_map <- c(
  C = "Control (C)",
  K = "Xan (K)"
)

label_group <- function(group) {
  group <- as.character(group)
  out <- group
  hit <- group %in% names(group_label_map)
  out[hit] <- unname(group_label_map[group[hit]])
  out
}

label_contrast <- function(cname, numerator = "K", denominator = "C") {
  sprintf("%s = %s vs %s", cname, label_group(numerator), label_group(denominator))
}

safe_read <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    msg <- sprintf("[missing] %s", path)
    if (required) stop(msg) else {
      cat("[warn] ", msg, "\n", sep = "")
      return(data.table())
    }
  }
  fread(path)
}

nonzero_ratio <- function(num, den) {
  out <- rep(NA_real_, length(num))
  ok <- !is.na(num) & !is.na(den) & den != 0
  out[ok] <- num[ok] / den[ok]
  out
}

classify_delta <- function(delta) {
  fifelse(is.na(delta), "not_comparable",
          fifelse(delta > 0, "higher_in_xan",
                  fifelse(delta < 0, "higher_in_control", "no_difference")))
}

module_root <- here::here()
project_root <- normalizePath(file.path(module_root, "..", ".."))
config_path <- file.path(module_root, "configs", "cellchat_params.yaml")
if (!file.exists(config_path)) stop("[CONTRACT] Config not found: ", config_path)
cfg <- yaml::read_yaml(config_path)

dataset_id <- Sys.getenv("DATASET_ID", "")
if (dataset_id == "" && !is.null(cfg$dataset$id)) dataset_id <- cfg$dataset$id
if (dataset_id == "") stop("[CONTRACT] No dataset_id: set DATASET_ID env var or config dataset.id")

base_resolved <- gsub("\\{dataset_id\\}", dataset_id, cfg$output$base_dir)
out_base <- file.path(project_root, base_resolved)

paths <- list(
  base = out_base,
  qc = file.path(out_base, "qc"),
  reports = file.path(out_base, "reports"),
  comparison_reports = file.path(out_base, "comparison", "reports"),
  deep_reports = file.path(out_base, "deep_dive", "reports"),
  summary_reports = file.path(out_base, "summary", "reports")
)
dir.create(paths$summary_reports, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("[dataset] %s\n", dataset_id))
cat(sprintf("[paths] out_base: %s\n", out_base))

contrast <- "K_vs_C"
numerator <- "K"
denominator <- "C"
contrast_label <- label_contrast(contrast, numerator, denominator)
groups <- c("C", "K")

written <- list()
write_report <- function(dt, filename, sources) {
  out_path <- file.path(paths$summary_reports, filename)
  fwrite(dt, out_path)
  cat(sprintf("[write] %s (%d rows)\n", out_path, nrow(dt)))
  written[[length(written) + 1L]] <<- data.table(
    file = filename,
    rows = nrow(dt),
    sources = paste(sources, collapse = ";"),
    generated_at = as.character(Sys.time()),
    dataset_id = dataset_id
  )
}

# 01 overview ---------------------------------------------------------------
inference <- safe_read(file.path(paths$qc, "inference_summary.csv"))
inference[, group_label := label_group(group)]
cell_counts <- safe_read(file.path(paths$qc, "celltype_group_counts.csv"), required = FALSE)

lr_rows <- rbindlist(lapply(groups, function(g) {
  dt <- safe_read(file.path(paths$reports, sprintf("LR_communication_%s.csv", g)), required = FALSE)
  data.table(group = g, lr_report_rows = nrow(dt))
}), fill = TRUE)
pathway_rows <- rbindlist(lapply(groups, function(g) {
  dt <- safe_read(file.path(paths$reports, sprintf("pathway_communication_%s.csv", g)), required = FALSE)
  data.table(group = g, pathway_report_rows = nrow(dt))
}), fill = TRUE)
cell_count_long <- data.table(group = groups)
if (nrow(cell_counts) > 0) {
  cell_count_long <- melt(cell_counts, id.vars = "celltype", variable.name = "group", value.name = "n_cells_by_celltype")
  cell_count_long <- cell_count_long[, .(
    celltypes_in_counts = .N,
    total_cells_in_counts = sum(n_cells_by_celltype, na.rm = TRUE)
  ), by = group]
}
overview <- merge(inference, lr_rows, by = "group", all.x = TRUE)
overview <- merge(overview, pathway_rows, by = "group", all.x = TRUE)
overview <- merge(overview, cell_count_long, by = "group", all.x = TRUE)
overview[, dataset_id := dataset_id]
setcolorder(overview, c("dataset_id", "group", "group_label"))
write_report(
  overview,
  "01_run_overview.csv",
  c("qc/inference_summary.csv", "qc/celltype_group_counts.csv", "reports/LR_communication_{C,K}.csv", "reports/pathway_communication_{C,K}.csv")
)

# 02 pathway presence -------------------------------------------------------
presence <- safe_read(file.path(paths$comparison_reports, "pathway_presence_matrix.csv"))
presence[, `:=`(
  dataset_id = dataset_id,
  contrast = contrast,
  contrast_label = contrast_label,
  present_control = as.integer(C %in% 1),
  present_xan = as.integer(K %in% 1)
)]
presence[, presence_category := fifelse(C == 1 & K == 1, "shared",
                                        fifelse(C == 0 & K == 1, "Xan (K)-only",
                                                fifelse(C == 1 & K == 0, "Control (C)-only", "absent")))]
setcolorder(presence, c("dataset_id", "contrast", "contrast_label", "pathway", "presence_category", "present_control", "present_xan", "C", "K", "n_groups"))
write_report(presence, "02_pathway_presence_summary.csv", "comparison/reports/pathway_presence_matrix.csv")

# 03 pathway communication --------------------------------------------------
pathway_long <- rbindlist(lapply(groups, function(g) {
  dt <- safe_read(file.path(paths$reports, sprintf("pathway_communication_%s.csv", g)), required = FALSE)
  if (nrow(dt) == 0) return(data.table())
  agg <- dt[, .(
    n_source_target_edges = .N,
    total_prob = sum(prob, na.rm = TRUE),
    mean_prob = mean(prob, na.rm = TRUE),
    max_prob = max(prob, na.rm = TRUE),
    min_pval = min(pval, na.rm = TRUE)
  ), by = .(pathway = pathway_name)]
  agg[, group := g]
  agg[, group_label := label_group(g)]
  agg
}), fill = TRUE)

if (nrow(pathway_long) > 0) {
  pathway_wide <- dcast(pathway_long, pathway ~ group, value.var = "total_prob")
  if (!"C" %in% names(pathway_wide)) pathway_wide[, C := NA_real_]
  if (!"K" %in% names(pathway_wide)) pathway_wide[, K := NA_real_]
  setnames(pathway_wide, c("C", "K"), c("total_prob_control", "total_prob_xan"))
  pathway_wide[, `:=`(
    delta_total_prob_K_minus_C = total_prob_xan - total_prob_control,
    ratio_total_prob_K_over_C = nonzero_ratio(total_prob_xan, total_prob_control),
    dominant_direction = classify_delta(total_prob_xan - total_prob_control),
    contrast = contrast,
    contrast_label = contrast_label,
    dataset_id = dataset_id
  )]
  pathway_edges <- pathway_long[, .(
    n_edges_control = sum(fifelse(group == "C", n_source_target_edges, 0L), na.rm = TRUE),
    n_edges_xan = sum(fifelse(group == "K", n_source_target_edges, 0L), na.rm = TRUE)
  ), by = pathway]
  pathway_summary <- merge(pathway_wide, pathway_edges, by = "pathway", all.x = TRUE)
  pathway_summary <- merge(pathway_summary, presence[, .(pathway, presence_category)], by = "pathway", all.x = TRUE)
  setcolorder(pathway_summary, c("dataset_id", "contrast", "contrast_label", "pathway", "presence_category"))
} else {
  pathway_summary <- data.table()
}
write_report(pathway_summary, "03_pathway_communication_summary.csv", "reports/pathway_communication_{C,K}.csv")

# 04 LR communication -------------------------------------------------------
lr_summary <- rbindlist(lapply(groups, function(g) {
  dt <- safe_read(file.path(paths$reports, sprintf("LR_communication_%s.csv", g)), required = FALSE)
  if (nrow(dt) == 0) return(data.table())
  annotation_col <- if ("annotation" %in% names(dt)) "annotation" else NULL
  if (is.null(annotation_col)) dt[, annotation := NA_character_]
  interaction_col <- if ("interaction_name" %in% names(dt)) "interaction_name" else NULL
  dt[, source_target := paste(source, target, sep = " -> ")]
  dt[, interaction_id := if (!is.null(interaction_col)) get(interaction_col) else NA_character_]
  agg <- dt[, .(
    n_lr_rows = .N,
    n_unique_interactions = uniqueN(interaction_id, na.rm = TRUE),
    n_source_target_pairs = uniqueN(source_target),
    total_prob = sum(prob, na.rm = TRUE),
    mean_prob = mean(prob, na.rm = TRUE),
    max_prob = max(prob, na.rm = TRUE),
    min_pval = min(pval, na.rm = TRUE)
  ), by = .(pathway = pathway_name, annotation)]
  agg[, group := g]
  agg[, group_label := label_group(g)]
  agg
}), fill = TRUE)
lr_summary[, `:=`(dataset_id = dataset_id, contrast = contrast, contrast_label = contrast_label)]
setcolorder(lr_summary, c("dataset_id", "contrast", "contrast_label", "group", "group_label"))
write_report(lr_summary, "04_lr_communication_summary.csv", "reports/LR_communication_{C,K}.csv")

# 05 celltype pair summary --------------------------------------------------
pair_tables <- lapply(groups, function(g) {
  dt <- safe_read(file.path(paths$comparison_reports, sprintf("celltype_pairs_%s.csv", g)), required = FALSE)
  if (nrow(dt) == 0) return(data.table(sender = character(), receiver = character()))
  setnames(dt, c("n_interactions", "interaction_weight"), sprintf(c("n_interactions_%s", "interaction_weight_%s"), g))
  dt
})
pair_summary <- Reduce(function(x, y) merge(x, y, by = c("sender", "receiver"), all = TRUE), pair_tables)
if (nrow(pair_summary) > 0) {
  for (nm in c("n_interactions_C", "n_interactions_K", "interaction_weight_C", "interaction_weight_K")) {
    if (!nm %in% names(pair_summary)) pair_summary[, (nm) := 0]
    pair_summary[is.na(get(nm)), (nm) := 0]
  }
  pair_summary[, `:=`(
    delta_n_interactions_K_minus_C = n_interactions_K - n_interactions_C,
    delta_interaction_weight_K_minus_C = interaction_weight_K - interaction_weight_C,
    ratio_interaction_weight_K_over_C = nonzero_ratio(interaction_weight_K, interaction_weight_C),
    dominant_direction = classify_delta(interaction_weight_K - interaction_weight_C),
    dataset_id = dataset_id,
    contrast = contrast,
    contrast_label = contrast_label
  )]
  pair_summary[, abs_delta_weight := abs(delta_interaction_weight_K_minus_C)]
  setorder(pair_summary, -abs_delta_weight)
  pair_summary[, abs_delta_weight := NULL]
  setcolorder(pair_summary, c("dataset_id", "contrast", "contrast_label", "sender", "receiver"))
}
write_report(pair_summary, "05_celltype_pair_summary.csv", "comparison/reports/celltype_pairs_{C,K}.csv")

# 06 focused pair probability -----------------------------------------------
pair_prob_files <- list.files(paths$deep_reports, pattern = "^07_pair_prob_.*\\.csv$", full.names = TRUE)
focused_pairs <- rbindlist(lapply(pair_prob_files, function(path) {
  dt <- safe_read(path, required = FALSE)
  if (nrow(dt) == 0) return(data.table())
  dt[, pair_file := basename(path)]
  if (!"group_label" %in% names(dt)) dt[, group_label := label_group(group)]
  dt
}), fill = TRUE)
if (nrow(focused_pairs) > 0) {
  focused_wide <- dcast(focused_pairs, pair_file + source + target ~ group, value.var = "total_prob", fun.aggregate = sum)
  if (!"C" %in% names(focused_wide)) focused_wide[, C := NA_real_]
  if (!"K" %in% names(focused_wide)) focused_wide[, K := NA_real_]
  setnames(focused_wide, c("C", "K"), c("total_prob_control", "total_prob_xan"))
  focused_wide[, `:=`(
    delta_total_prob_K_minus_C = total_prob_xan - total_prob_control,
    ratio_total_prob_K_over_C = nonzero_ratio(total_prob_xan, total_prob_control),
    dominant_direction = classify_delta(total_prob_xan - total_prob_control),
    dataset_id = dataset_id,
    contrast = contrast,
    contrast_label = contrast_label
  )]
  setorder(focused_wide, delta_total_prob_K_minus_C)
  setcolorder(focused_wide, c("dataset_id", "contrast", "contrast_label", "pair_file", "source", "target"))
} else {
  focused_wide <- data.table()
}
write_report(focused_wide, "06_deep_dive_pair_probability_summary.csv", "deep_dive/reports/07_pair_prob_*.csv")

# 07 LR contribution --------------------------------------------------------
lr_contrib_files <- list.files(paths$deep_reports, pattern = "^03_lr_contribution_.*\\.csv$", full.names = TRUE)
lr_contrib <- rbindlist(lapply(lr_contrib_files, function(path) {
  dt <- safe_read(path, required = FALSE)
  if (nrow(dt) == 0) return(data.table())
  dt[, contribution_file := basename(path)]
  if (!"group_label" %in% names(dt)) dt[, group_label := label_group(group)]
  dt
}), fill = TRUE)
if (nrow(lr_contrib) > 0) {
  lr_contrib_summary <- lr_contrib[, .(
    contribution = sum(contribution, na.rm = TRUE),
    mean_contribution = mean(contribution, na.rm = TRUE),
    n_rows = .N
  ), by = .(contribution_file, pathway, name, group, group_label)]
  lr_contrib_summary[, `:=`(dataset_id = dataset_id, contrast = contrast, contrast_label = contrast_label)]
  setcolorder(lr_contrib_summary, c("dataset_id", "contrast", "contrast_label", "contribution_file", "pathway", "name", "group", "group_label"))
} else {
  lr_contrib_summary <- data.table()
}
write_report(lr_contrib_summary, "07_lr_contribution_summary.csv", "deep_dive/reports/03_lr_contribution_*.csv")

index <- rbindlist(written, fill = TRUE)
write_report(index, "00_summary_index.csv", "summary script generated files")

cat("\n")
cat("+================================================================+\n")
cat("| CellChat Step 05 Complete\n")
cat("+================================================================+\n")
cat(sprintf("| Dataset: %s\n", dataset_id))
cat(sprintf("| Summary: %s\n", paths$summary_reports))
cat(sprintf("| Files:   %d\n", nrow(index)))
cat("+================================================================+\n")
