# ============================================================================
# 03_compare_groups.R -- CellChat group comparison (config-driven)
#
# Purpose:  Merge per-group inferred CellChat objects and run pairwise/global
#           comparison plots and summary exports.
# Input:    results/scrna/{dataset_id}/09_cellchat/objects/
# Output:   results/scrna/{dataset_id}/09_cellchat/comparison/{plots,reports}
# Usage:    cd modules/cellchat && Rscript scripts/03_compare_groups.R
# Override: DATASET_ID=gse253718_luad Rscript scripts/03_compare_groups.R
# Tested:   CellChat v2.x API
# ============================================================================

cat(
  "\n+================================================================+\n",
  "|               CellChat Step 03: Group Comparison              |\n",
  "+================================================================+\n",
  sep = ""
)

suppressPackageStartupMessages({
  library(CellChat)
  library(patchwork)
  library(ggplot2)
  library(yaml)
  library(here)
  library(data.table)
})

# ══════════════════════════════════════════════════════════════════════════════
# Preflight: verify critical packages are loadable
# ══════════════════════════════════════════════════════════════════════════════
required_pkgs <- c("CellChat", "Seurat", "yaml", "data.table", "patchwork", "ggplot2")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(sprintf(
    "[PREFLIGHT] Missing packages: %s\n  → Fix: cd modules/cellchat && Rscript --vanilla -e 'renv::hydrate(); renv::snapshot()'",
    paste(missing_pkgs, collapse = ", ")
  ))
}



t0 <- Sys.time()

# ★ FIX 5: Global plot failure counter (accumulated across all sections)
plot_section_total  <- 0L
plot_section_failed <- 0L

# Provide `%||%` if not available in the current runtime.
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) {
    if (is.null(a)) b else a
  }
}

safe_dev_off <- function() {
  try(dev.off(), silent = TRUE)
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

label_contrast <- function(cname, numerator = NULL, denominator = NULL) {
  if (!is.null(numerator) && !is.null(denominator)) {
    return(sprintf("%s = %s vs %s", cname, label_group(numerator), label_group(denominator)))
  }
  cname
}

normalize_comparisons <- function(comparisons) {
  comp_names <- names(comparisons)
  out <- vector("list", length(comparisons))

  for (i in seq_along(comparisons)) {
    comp <- comparisons[[i]]
    fallback_name <- if (!is.null(comp_names) && nzchar(comp_names[i])) {
      comp_names[i]
    } else {
      sprintf("comparison_%02d", i)
    }

    comp_name <- as.character(comp$name %||% fallback_name)

    if (!is.null(comp$groups) && length(comp$groups) >= 2) {
      numerator <- as.character(comp$groups[1])
      denominator <- as.character(comp$groups[2])
    } else {
      numerator <- as.character(comp$numerator %||% NA_character_)
      denominator <- as.character(comp$denominator %||% NA_character_)
    }

    if (is.na(numerator) || is.na(denominator) || numerator == "" || denominator == "") {
      stop(sprintf(
        "[CONTRACT] Invalid comparison '%s': need groups[1:2] or numerator/denominator.",
        comp_name
      ))
    }

    out[[i]] <- list(
      name = comp_name,
      numerator = numerator,
      denominator = denominator
    )
  }

  out
}

load_cellchat_objects <- function(obj_dir, analysis_groups) {
  list_inferred <- file.path(obj_dir, "cellchat_list_inferred.rds")
  if (file.exists(list_inferred)) {
    cc_list <- readRDS(list_inferred)
    if (!is.list(cc_list)) {
      stop("[CONTRACT] cellchat_list_inferred.rds is not a list.")
    }

    missing_from_list <- setdiff(analysis_groups, names(cc_list))
    if (length(missing_from_list) == 0) {
      return(list(
        objects = cc_list[analysis_groups],
        source_mode = "cellchat_list_inferred"
      ))
    }

    cat(sprintf(
      "[WARN] cellchat_list_inferred.rds missing groups: %s; trying per-group fallbacks.\n",
      paste(missing_from_list, collapse = ", ")
    ))
  }

  per_group_inferred <- file.path(obj_dir, sprintf("cellchat_%s_inferred.rds", analysis_groups))
  if (all(file.exists(per_group_inferred))) {
    cc_list <- lapply(per_group_inferred, readRDS)
    names(cc_list) <- analysis_groups
    return(list(
      objects = cc_list,
      source_mode = "per_group_inferred"
    ))
  }

  per_group_legacy <- file.path(obj_dir, sprintf("cellchat_%s.rds", analysis_groups))
  if (all(file.exists(per_group_legacy))) {
    cc_list <- lapply(per_group_legacy, readRDS)
    names(cc_list) <- analysis_groups
    return(list(
      objects = cc_list,
      source_mode = "per_group_legacy"
    ))
  }

  missing_inferred <- analysis_groups[!file.exists(per_group_inferred)]
  missing_legacy <- analysis_groups[!file.exists(per_group_legacy)]

  stop(sprintf(
    paste0(
      "[CONTRACT] Cannot resolve Step 03 input objects for all analysis groups.\n",
      "  Tried 1) cellchat_list_inferred.rds\n",
      "  Tried 2) cellchat_{group}_inferred.rds (missing: %s)\n",
      "  Tried 3) cellchat_{group}.rds (missing: %s)"
    ),
    if (length(missing_inferred) == 0) "none" else paste(missing_inferred, collapse = ", "),
    if (length(missing_legacy) == 0) "none" else paste(missing_legacy, collapse = ", ")
  ))
}

# ============================================================================
# 0. Bootstrap config and paths
# ============================================================================
module_root <- here::here()
project_root <- normalizePath(file.path(module_root, "..", ".."))

config_path <- file.path(module_root, "configs", "cellchat_params.yaml")
if (!file.exists(config_path)) {
  stop("[CONTRACT] Config not found: ", config_path)
}
cfg <- yaml::read_yaml(config_path)

dataset_id <- Sys.getenv("DATASET_ID", "")
if (dataset_id == "" && !is.null(cfg$dataset$id)) {
  dataset_id <- cfg$dataset$id
}
if (dataset_id == "") {
  stop("[CONTRACT] No dataset_id: set DATASET_ID env var or config dataset.id")
}

base_template <- cfg$output$base_dir
base_resolved <- gsub("\\{dataset_id\\}", dataset_id, base_template)
out_base <- file.path(project_root, base_resolved)

out <- list(
  base = out_base,
  objects = file.path(out_base, "objects"),
  qc = file.path(out_base, "qc"),
  comparison = file.path(out_base, "comparison"),
  plots = file.path(out_base, "comparison", "plots"),
  reports = file.path(out_base, "comparison", "reports")
)

for (d in out[c("objects", "qc", "comparison", "plots", "reports")]) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat(sprintf("[config] Loaded: %s\n", config_path))
cat(sprintf("[dataset] id: %s\n", dataset_id))
cat(sprintf("[paths] module_root:  %s\n", module_root))
cat(sprintf("[paths] project_root: %s\n", project_root))
cat(sprintf("[paths] out_base:     %s\n", out$base))
cat(sprintf("[env] CellChat %s\n", as.character(packageVersion("CellChat"))))

# Comparison-specific optional config with defaults.
default_focus_pathways <- c(
  "SPP1", "MIF", "CCL", "CXCL", "CSF3", "GALECTIN",
  "TNF", "IFN-II", "IL6", "IL10", "CD86", "ICAM",
  "MHC-I", "MHC-II", "PD-L1", "CD80", "COMPLEMENT",
  "FN1", "COLLAGEN", "LAMININ", "VEGF", "PDGF",
  "WNT", "NOTCH", "TGFb", "BMP", "EGF", "GAS"
)

comparison_cfg <- cfg$comparison %||% list()

focus_pathways <- comparison_cfg$focus_pathways %||% default_focus_pathways
focus_source <- if (is.null(comparison_cfg$focus_pathways)) "default" else "config"
focus_pathways <- unique(as.character(focus_pathways))
if (length(focus_pathways) == 0) {
  focus_pathways <- default_focus_pathways
  focus_source <- "default"
}

bubble_angle_x <- comparison_cfg$bubble_angle_x %||% 45
if (!is.numeric(bubble_angle_x) || length(bubble_angle_x) != 1 || is.na(bubble_angle_x)) {
  bubble_angle_x <- 45
}

bubble_remove_isolate <- comparison_cfg$bubble_remove_isolate
if (is.null(bubble_remove_isolate)) {
  bubble_remove_isolate <- TRUE
} else {
  bubble_remove_isolate <- isTRUE(bubble_remove_isolate)
}

cat(sprintf("[config] comparison.focus_pathways: %d (%s)\n", length(focus_pathways), focus_source))
cat(sprintf("[config] comparison.bubble_angle_x: %s\n", bubble_angle_x))
cat(sprintf("[config] comparison.bubble_remove_isolate: %s\n", bubble_remove_isolate))

# ============================================================================
# 1. Load design and normalize analysis scope
# ============================================================================
design_path <- file.path(out$objects, "design.rds")
if (!file.exists(design_path)) {
  stop("[CONTRACT] design.rds not found. Run 01_prepare_cellchat.R first.\n  Expected: ", design_path)
}
design <- readRDS(design_path)

if (is.null(design$all_groups) || length(design$all_groups) == 0) {
  stop("[CONTRACT] design.rds missing non-empty all_groups")
}

analysis_groups <- unique(as.character(design$all_groups))
reference_group <- design$reference_group %||% NULL
if (!is.null(reference_group) && reference_group %in% analysis_groups) {
  analysis_groups <- c(reference_group, setdiff(analysis_groups, reference_group))
}

if (is.null(design$comparisons) || length(design$comparisons) == 0) {
  stop("[CONTRACT] design.rds missing non-empty comparisons")
}

comparisons <- normalize_comparisons(design$comparisons)

bad_comparisons <- character(0)
for (comp in comparisons) {
  if (!(comp$numerator %in% analysis_groups) || !(comp$denominator %in% analysis_groups)) {
    bad_comparisons <- c(
      bad_comparisons,
      sprintf("%s(%s vs %s)", comp$name, comp$numerator, comp$denominator)
    )
  }
}
if (length(bad_comparisons) > 0) {
  stop(sprintf(
    "[CONTRACT] Comparisons reference groups outside analysis scope: %s\n  analysis_groups: %s",
    paste(bad_comparisons, collapse = ", "),
    paste(analysis_groups, collapse = ", ")
  ))
}

cat(sprintf("[load] Design groups (%d): %s\n", length(analysis_groups), paste(analysis_groups, collapse = ", ")))
cat(sprintf("[load] Comparisons (%d): %s\n",
            length(comparisons), paste(vapply(comparisons, `[[`, "", "name"), collapse = ", ")))

# ============================================================================
# 2. Load CellChat objects using fallback chain
# ============================================================================
loaded <- load_cellchat_objects(out$objects, analysis_groups)
cc_list <- loaded$objects
source_mode <- loaded$source_mode
cat(sprintf("[load] Input mode: %s\n", source_mode))

for (grp in analysis_groups) {
  cc <- cc_list[[grp]]
  if (is.null(cc@DB) || length(cc@DB) == 0) {
    stop(sprintf("[CONTRACT] %s object missing @DB. Step 01/02 artifact invalid.", grp))
  }
  if (is.null(cc@idents) || length(cc@idents) == 0) {
    stop(sprintf("[CONTRACT] %s object missing @idents.", grp))
  }
  if (is.null(cc@net) || length(cc@net) == 0 || is.null(cc@netP) || length(cc@netP) == 0) {
    stop(sprintf(
      "[CONTRACT] %s object lacks inferred slots (@net/@netP). Run 02_run_cellchat.R first.",
      grp
    ))
  }

  n_pathways <- length(cc@netP$pathways %||% character(0))
  cat(sprintf("[load] %s: %d cells, %d celltypes, %d pathways\n",
              grp, length(cc@idents), length(levels(cc@idents)), n_pathways))
}

# ============================================================================
# 3. Merge groups
# ============================================================================
cat("\n[step] Merging groups...\n")

cc_merged <- tryCatch(
  mergeCellChat(cc_list, add.names = analysis_groups),
  error = function(e) {
    stop(sprintf(
      "[FATAL] mergeCellChat failed under CellChat %s: %s\n  Hint: API/signature may differ across versions.",
      as.character(packageVersion("CellChat")),
      e$message
    ))
  }
)

# ★ FIX 3: Validate merged dataset membership using metadata labels.
if (is.null(cc_merged@meta) || !"datasets" %in% colnames(cc_merged@meta)) {
  stop("[CONTRACT] Merged CellChat object missing meta$datasets label column")
}

merged_groups <- unique(as.character(cc_merged@meta$datasets))
n_merged <- length(merged_groups)

missing_in_merged <- setdiff(analysis_groups, merged_groups)
extra_in_merged <- setdiff(merged_groups, analysis_groups)

if (length(missing_in_merged) > 0 || length(extra_in_merged) > 0) {
  stop(sprintf(
    paste0(
      "[CONTRACT] Merged CellChat group mismatch.\n",
      "  expected: %s\n",
      "  observed: %s\n",
      "  missing:  %s\n",
      "  extra:    %s"
    ),
    paste(analysis_groups, collapse = ", "),
    paste(merged_groups, collapse = ", "),
    if (length(missing_in_merged) == 0) "none" else paste(missing_in_merged, collapse = ", "),
    if (length(extra_in_merged) == 0) "none" else paste(extra_in_merged, collapse = ", ")
  ))
}

cat(sprintf("[validate] Merged object contains %d datasets: %s\n",
            n_merged, paste(merged_groups, collapse = ", ")))

merged_dynamic <- file.path(out$objects, sprintf("cellchat_merged_%dgrp.rds", length(analysis_groups)))
merged_alias <- file.path(out$objects, "cellchat_merged.rds")
saveRDS(cc_merged, merged_dynamic)
saveRDS(cc_merged, merged_alias)
cat(sprintf("[save] %s\n", merged_dynamic))
cat(sprintf("[save] %s\n", merged_alias))

# ============================================================================
# 4. Global comparison: interaction count and strength
# ============================================================================
cat("\n[step] Global interaction comparison...\n")

plot_section_total <- plot_section_total + 1L

# ★ FIX 2: compareInteractions expects group names (character), not integer indices
tryCatch({
  p1 <- compareInteractions(cc_merged, show.legend = TRUE, group = analysis_groups) +
    labs(title = "Interaction count", subtitle = paste(label_group(analysis_groups), collapse = " vs "))
  p2 <- compareInteractions(cc_merged, show.legend = TRUE, group = analysis_groups, measure = "weight") +
    labs(title = "Interaction strength", subtitle = paste(label_group(analysis_groups), collapse = " vs "))

  pdf(file.path(out$plots, "01_global_interaction_comparison.pdf"), width = 12, height = 6)
  print(p1 + p2)
  safe_dev_off()

  png(file.path(out$plots, "01_global_interaction_comparison.png"), width = 1200, height = 600, res = 150)
  print(p1 + p2)
  safe_dev_off()

  cat("[plot] 01_global_interaction_comparison\n")
}, error = function(e) {
  cat(sprintf("[WARN] Global interaction comparison failed: %s\n", e$message))
  safe_dev_off()
  plot_section_failed <<- plot_section_failed + 1L
})

# ============================================================================
# 5. Pairwise differential interaction plots
# ============================================================================
cat("\n[step] Differential interaction plots (pairwise)...\n")

for (comp in comparisons) {
  cname <- comp$name
  num <- comp$numerator
  den <- comp$denominator

  idx_num <- match(num, analysis_groups)
  idx_den <- match(den, analysis_groups)

  if (is.na(idx_num) || is.na(idx_den)) {
    cat(sprintf("[skip] %s: group not found in merged object (%s vs %s)\n", cname, num, den))
    next
  }

  contrast_label <- label_contrast(cname, num, den)
  cat(sprintf("\n  --- %s ---\n", contrast_label))

  plot_section_total <- plot_section_total + 1L
  tryCatch({
    pdf(file.path(out$plots, sprintf("02_diff_interaction_%s.pdf", cname)), width = 14, height = 6)
    par(mfrow = c(1, 2))
    netVisual_diffInteraction(
      cc_merged,
      comparison = c(idx_den, idx_num),
      weight.scale = TRUE,
      measure = "count",
      title.name = sprintf("Diff #interactions: %s", contrast_label)
    )
    netVisual_diffInteraction(
      cc_merged,
      comparison = c(idx_den, idx_num),
      weight.scale = TRUE,
      measure = "weight",
      title.name = sprintf("Diff strength: %s", contrast_label)
    )
    safe_dev_off()

    png(file.path(out$plots, sprintf("02_diff_interaction_%s.png", cname)), width = 1400, height = 600, res = 150)
    par(mfrow = c(1, 2))
    netVisual_diffInteraction(
      cc_merged,
      comparison = c(idx_den, idx_num),
      weight.scale = TRUE,
      measure = "count",
      title.name = sprintf("Diff #interactions: %s", contrast_label)
    )
    netVisual_diffInteraction(
      cc_merged,
      comparison = c(idx_den, idx_num),
      weight.scale = TRUE,
      measure = "weight",
      title.name = sprintf("Diff strength: %s", contrast_label)
    )
    safe_dev_off()

    cat(sprintf("  [plot] 02_diff_interaction_%s\n", cname))
  }, error = function(e) {
    cat(sprintf("  [WARN] Differential interaction plot failed for %s: %s\n", cname, e$message))
    safe_dev_off()
    plot_section_failed <<- plot_section_failed + 1L
  })

  plot_section_total <- plot_section_total + 1L
  tryCatch({
    pdf(file.path(out$plots, sprintf("03_diff_heatmap_%s.pdf", cname)), width = 12, height = 10)
    par(mfrow = c(1, 2))
    netVisual_heatmap(
      cc_merged,
      comparison = c(idx_den, idx_num),
      measure = "count",
      title.name = sprintf("#interactions: %s", contrast_label)
    )
    netVisual_heatmap(
      cc_merged,
      comparison = c(idx_den, idx_num),
      measure = "weight",
      title.name = sprintf("Strength: %s", contrast_label)
    )
    safe_dev_off()

    png(file.path(out$plots, sprintf("03_diff_heatmap_%s.png", cname)), width = 1200, height = 1000, res = 150)
    par(mfrow = c(1, 2))
    netVisual_heatmap(
      cc_merged,
      comparison = c(idx_den, idx_num),
      measure = "count",
      title.name = sprintf("#interactions: %s", contrast_label)
    )
    netVisual_heatmap(
      cc_merged,
      comparison = c(idx_den, idx_num),
      measure = "weight",
      title.name = sprintf("Strength: %s", contrast_label)
    )
    safe_dev_off()

    cat(sprintf("  [plot] 03_diff_heatmap_%s\n", cname))
  }, error = function(e) {
    cat(sprintf("  [WARN] Differential heatmap failed for %s: %s\n", cname, e$message))
    safe_dev_off()
    plot_section_failed <<- plot_section_failed + 1L
  })
}

# ============================================================================
# 6. Information flow comparison
# ============================================================================
cat("\n[step] Information flow comparison (rankNet)...\n")

plot_section_total <- plot_section_total + 1L
tryCatch({
  p_flow <- rankNet(cc_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE) +
    labs(title = sprintf("Information flow: %s", label_contrast("K_vs_C", "K", "C")))

  pdf(file.path(out$plots, "04_information_flow_comparison.pdf"), width = 10, height = 14)
  print(p_flow)
  safe_dev_off()

  png(file.path(out$plots, "04_information_flow_comparison.png"), width = 1000, height = 1400, res = 150)
  print(p_flow)
  safe_dev_off()

  cat("[plot] 04_information_flow_comparison\n")
}, error = function(e) {
  cat(sprintf("[WARN] rankNet stacked comparison failed: %s\n", e$message))
  safe_dev_off()
  plot_section_failed <<- plot_section_failed + 1L
})

plot_section_total <- plot_section_total + 1L
tryCatch({
  p_flow2 <- rankNet(cc_merged, mode = "comparison", stacked = FALSE, do.stat = TRUE) +
    labs(title = sprintf("Information flow (unstacked): %s", label_contrast("K_vs_C", "K", "C")))

  pdf(file.path(out$plots, "04b_information_flow_unstacked.pdf"), width = 10, height = 14)
  print(p_flow2)
  safe_dev_off()

  png(file.path(out$plots, "04b_information_flow_unstacked.png"), width = 1000, height = 1400, res = 150)
  print(p_flow2)
  safe_dev_off()

  cat("[plot] 04b_information_flow_unstacked\n")
}, error = function(e) {
  cat(sprintf("[WARN] rankNet unstacked comparison failed: %s\n", e$message))
  safe_dev_off()
  plot_section_failed <<- plot_section_failed + 1L
})

# ============================================================================
# 7. Pathway-level bubble plots (chunked to avoid overplotting)
# ============================================================================
cat("\n[step] Pathway-level bubble plots (chunked)...\n")

# ── Read bubble config ──
bubble_chunk_size <- comparison_cfg$bubble_chunk_size
if (is.null(bubble_chunk_size) || !is.numeric(bubble_chunk_size) || bubble_chunk_size < 1) {
  bubble_chunk_size <- 10L  # default: 10 pathways per page
}
bubble_chunk_size <- as.integer(bubble_chunk_size)

bubble_font_size <- comparison_cfg$bubble_font_size
if (is.null(bubble_font_size) || !is.numeric(bubble_font_size)) {
  bubble_font_size <- 10
}

bubble_thresh <- comparison_cfg$bubble_thresh  # NULL = show all
if (!is.null(bubble_thresh) && (!is.numeric(bubble_thresh) || is.na(bubble_thresh))) {
  bubble_thresh <- NULL
}

cat(sprintf("[config] bubble: chunk_size=%d, font_size=%s, thresh=%s\n",
            bubble_chunk_size, bubble_font_size,
            if (is.null(bubble_thresh)) "NULL (all)" else as.character(bubble_thresh)))

# ── Helper: split vector into chunks ──
chunk_vector <- function(x, size) {
  split(x, ceiling(seq_along(x) / size))
}

# ============================================================================
# Adaptive bubble wrapper: auto-scale font + figure size based on axis density
# ============================================================================
adaptive_bubble <- function(cc_merged, comparison, signaling = NULL,
                            angle.x = 45, remove.isolate = TRUE,
                            thresh = NULL, title.name = "",
                            font.size = 10,
                            min_font = 4, max_font = 12,
                            px_per_x = 30, px_per_y = 18,
                            min_w = 8, max_w = 40,
                            min_h = 6, max_h = 60,
                            ...) {
  # ── 1. Generate ggplot object ──
  args <- list(
    object         = cc_merged,
    comparison     = comparison,
    angle.x        = angle.x,
    remove.isolate = remove.isolate,
    font.size      = font.size,
    title.name     = title.name
  )
  if (!is.null(signaling))  args$signaling <- signaling
  if (!is.null(thresh))     args$thresh    <- thresh

  p <- do.call(netVisual_bubble, c(args, list(...)))

  # ── 2. Probe actual axis element count from ggplot internals ──
  n_x <- n_y <- NA_integer_
  tryCatch({
    gb <- ggplot_build(p)
    # Discrete axes: count unique positions in first data layer
    if (length(gb$data) > 0 && nrow(gb$data[[1]]) > 0) {
      d1 <- gb$data[[1]]
      if ("x" %in% names(d1)) n_x <- length(unique(d1$x))
      if ("y" %in% names(d1)) n_y <- length(unique(d1$y))
    }
    # Fallback: panel_params tick labels
    if (is.na(n_x) || is.na(n_y)) {
      pp <- gb$layout$panel_params[[1]]
      if (is.na(n_x)) n_x <- length(pp$x$get_labels %||% pp$x.labels %||% character(0))
      if (is.na(n_y)) n_y <- length(pp$y$get_labels %||% pp$y.labels %||% character(0))
    }
  }, error = function(e) {
    # If introspection fails, use generous defaults
    n_x <<- 50L
    n_y <<- 50L
  })
  if (is.na(n_x)) n_x <- 50L
  if (is.na(n_y)) n_y <- 50L

  # ── 3. Dynamic font size: inversely proportional to axis density ──
  n_max_axis <- max(n_x, n_y, 1L)
  auto_font <- if (n_max_axis <= 20) {
    max_font
  } else if (n_max_axis <= 50) {
    round(max_font - (n_max_axis - 20) * (max_font - min_font) / 80, 1)
  } else {
    max(min_font, round(300 / n_max_axis, 1))
  }
  auto_font <- max(min_font, min(max_font, auto_font))

  # ── 4. Dynamic figure dimensions (inches) ──
  auto_w <- min(max_w, max(min_w, n_x * px_per_x / 100 + 4))
  auto_h <- min(max_h, max(min_h, n_y * px_per_y / 100 + 3))

  # ── 5. Dynamic x-label angle: go vertical if very crowded ──
  auto_angle <- if (n_x > 60) 90 else if (n_x > 30) 75 else angle.x
  auto_hjust <- if (auto_angle >= 75) 1 else if (auto_angle >= 30) 1 else 0.5
  auto_vjust <- if (auto_angle == 90) 0.5 else 1

  # ── 6. Apply theme overrides to the existing ggplot ──
  p <- p + theme(
    axis.text.x = element_text(
      size  = auto_font,
      angle = auto_angle,
      hjust = auto_hjust,
      vjust = auto_vjust
    ),
    axis.text.y = element_text(size = auto_font),
    legend.text = element_text(size = max(auto_font, 7)),
    plot.title  = element_text(size = auto_font + 2)
  )

  list(
    plot      = p,
    width_in  = auto_w,
    height_in = auto_h,
    n_x       = n_x,
    n_y       = n_y,
    font_size = auto_font,
    angle_x   = auto_angle
  )
}


for (comp in comparisons) {
  cname <- comp$name
  num   <- comp$numerator
  den   <- comp$denominator

  idx_num <- match(num, analysis_groups)
  idx_den <- match(den, analysis_groups)

  if (is.na(idx_num) || is.na(idx_den)) {
    cat(sprintf("  [skip] %s: group missing (%s vs %s)\n", cname, num, den))
    next
  }

  pw_num <- cc_list[[num]]@netP$pathways %||% character(0)
  pw_den <- cc_list[[den]]@netP$pathways %||% character(0)
  common_pw <- intersect(pw_num, pw_den)
  pw_plot   <- intersect(focus_pathways, common_pw)

  if (length(pw_plot) == 0) {
    cat(sprintf("  [skip] %s: no configured focus pathways in common\n", cname))
    next
  }

  cat(sprintf("  %s: %d/%d focus pathways present → %d chunk(s)\n",
              cname, length(pw_plot), length(focus_pathways),
              ceiling(length(pw_plot) / bubble_chunk_size)))



  # ── 7a. Full overview (all focus pathways, auto-scaled) ──
  plot_section_total <- plot_section_total + 1L
  tryCatch({
    ab <- adaptive_bubble(
      cc_merged,
      comparison     = c(idx_den, idx_num),
      signaling      = pw_plot,
      angle.x        = bubble_angle_x,
      remove.isolate = bubble_remove_isolate,
      thresh         = bubble_thresh,
      font.size      = bubble_font_size,
      title.name     = sprintf("All focus pathways: %s (%d pathways)", contrast_label, length(pw_plot))
    )

    cat(sprintf("  full: %d×%d axis elements, font=%.0f, angle=%d, %.0f×%.0f in\n",
                ab$n_x, ab$n_y, ab$font_size, ab$angle_x, ab$width_in, ab$height_in))

    pdf(file.path(out$plots, sprintf("05_bubble_focus_%s_full.pdf", cname)),
        width = ab$width_in, height = ab$height_in)
    print(ab$plot)
    safe_dev_off()

    png(file.path(out$plots, sprintf("05_bubble_focus_%s_full.png", cname)),
        width = ab$width_in * 100, height = ab$height_in * 100, res = 150)
    print(ab$plot)
    safe_dev_off()

    cat(sprintf("  [plot] 05_bubble_focus_%s_full\n", cname))
  }, error = function(e) {
    cat(sprintf("  [WARN] Full bubble plot failed for %s: %s\n", cname, e$message))
    safe_dev_off()
    plot_section_failed <<- plot_section_failed + 1L
  })


  # ── 7b. Chunked panels (readable per-page) ──
  chunks <- chunk_vector(pw_plot, bubble_chunk_size)

  for (ci in seq_along(chunks)) {
    chunk_pw <- chunks[[ci]]

    plot_section_total <- plot_section_total + 1L
    tryCatch({
      ab <- adaptive_bubble(
        cc_merged,
        comparison     = c(idx_den, idx_num),
        signaling      = chunk_pw,
        angle.x        = bubble_angle_x,
        remove.isolate = bubble_remove_isolate,
        thresh         = bubble_thresh,
        font.size      = bubble_font_size,
        title.name     = sprintf("%s chunk %d/%d (%d pathways)",
                                 contrast_label, ci, length(chunks), length(chunk_pw))
      )

      cat(sprintf("    chunk%02d: %d×%d axis elements, font=%.0f, angle=%d, %.0f×%.0f in\n",
                  ci, ab$n_x, ab$n_y, ab$font_size, ab$angle_x,
                  ab$width_in, ab$height_in))

      pdf(file.path(out$plots, sprintf("05_bubble_focus_%s_chunk%02d.pdf", cname, ci)),
          width = ab$width_in, height = ab$height_in)
      print(ab$plot)
      safe_dev_off()

      png(file.path(out$plots, sprintf("05_bubble_focus_%s_chunk%02d.png", cname, ci)),
          width  = ab$width_in * 100,
          height = ab$height_in * 100,
          res = 150)
      print(ab$plot)
      safe_dev_off()

      cat(sprintf("  [plot] 05_bubble_focus_%s_chunk%02d\n", cname, ci))
    }, error = function(e) {
      cat(sprintf("  [WARN] Bubble chunk %d failed for %s: %s\n", ci, cname, e$message))
      safe_dev_off()
      plot_section_failed <<- plot_section_failed + 1L
    })
  }


  # ── 7c. All common pathways (not just focus) — single overview ──
  if (length(common_pw) > length(pw_plot)) {
    plot_section_total <- plot_section_total + 1L
    tryCatch({
      ab_all <- adaptive_bubble(
        cc_merged,
        comparison     = c(idx_den, idx_num),
        angle.x        = bubble_angle_x,
        remove.isolate = bubble_remove_isolate,
        thresh         = bubble_thresh,
        font.size      = bubble_font_size,
        title.name     = sprintf("All common pathways: %s (%d pathways)", contrast_label, length(common_pw))
      )

      cat(sprintf("  all: %d×%d axis elements, font=%.0f, angle=%d, %.0f×%.0f in\n",
                  ab_all$n_x, ab_all$n_y, ab_all$font_size, ab_all$angle_x,
                  ab_all$width_in, ab_all$height_in))

      pdf(file.path(out$plots, sprintf("05_bubble_all_%s.pdf", cname)),
          width = ab_all$width_in, height = ab_all$height_in)
      print(ab_all$plot)
      safe_dev_off()

      png(file.path(out$plots, sprintf("05_bubble_all_%s.png", cname)),
          width  = ab_all$width_in * 100,
          height = ab_all$height_in * 100,
          res = 150)
      print(ab_all$plot)
      safe_dev_off()

      cat(sprintf("  [plot] 05_bubble_all_%s (%d pathways)\n", cname, length(common_pw)))
    }, error = function(e) {
      cat(sprintf("  [WARN] All-pathway bubble failed for %s: %s\n", cname, e$message))
      safe_dev_off()
      plot_section_failed <<- plot_section_failed + 1L
    })
  }

}


# ============================================================================
# 8. Signaling role comparison
# ============================================================================
cat("\n[step] Signaling role analysis (outgoing vs incoming)...\n")

# ★ FIX 4: Compute centrality on INDIVIDUAL objects BEFORE merging pairs.
#   netAnalysis_computeCentrality does not work on merged CellChat objects.
cat("  Computing signaling centrality per group...\n")
centrality_ok <- character(0)
for (grp in analysis_groups) {
  tryCatch({
    cc_list[[grp]] <- netAnalysis_computeCentrality(cc_list[[grp]])
    centrality_ok <- c(centrality_ok, grp)
    cat(sprintf("    [done] %s centrality computed\n", grp))
  }, error = function(e) {
    cat(sprintf("    [WARN] %s centrality FAILED: %s\n", grp, e$message))
  })
}
cat(sprintf("  Centrality computed for %d/%d groups\n",
            length(centrality_ok), length(analysis_groups)))

role_plot_total <- 0L
role_plot_failed <- 0L

for (comp in comparisons) {
  cname <- comp$name
  num <- comp$numerator
  den <- comp$denominator

  role_plot_total <- role_plot_total + 1L
  plot_section_total <- plot_section_total + 1L

  # Pre-check: both groups must have centrality computed
  if (!all(c(num, den) %in% centrality_ok)) {
    missing_cent <- setdiff(c(num, den), centrality_ok)
    cat(sprintf("  [SKIP] %s: centrality missing for %s — signaling role plot skipped\n",
                cname, paste(missing_cent, collapse = ", ")))
    plot_section_failed <- plot_section_failed + 1L
    role_plot_failed <- role_plot_failed + 1L
    next
  }

    tryCatch({
    contrast_label <- label_contrast(cname, num, den)
    cat(sprintf("  --- %s ---\n", contrast_label))

    # ★ FIX: CellChat v2 signaling role scatter is called on INDIVIDUAL objects,
    #   not merged. Each object must have centrality pre-computed (done above).
    p_den <- netAnalysis_signalingRole_scatter(
      cc_list[[den]],
      title = label_group(den)
    )
    p_num <- netAnalysis_signalingRole_scatter(
      cc_list[[num]],
      title = label_group(num)
    )
    p_combined <- p_den + p_num + plot_layout(ncol = 2)

    pdf(file.path(out$plots, sprintf("06_signaling_role_%s.pdf", cname)),
        width = 14, height = 7)
    print(p_combined)
    safe_dev_off()

    png(file.path(out$plots, sprintf("06_signaling_role_%s.png", cname)),
        width = 1400, height = 700, res = 150)
    print(p_combined)
    safe_dev_off()

    cat(sprintf("  [plot] 06_signaling_role_%s\n", cname))
  }, error = function(e) {
    cat(sprintf("  [WARN] Signaling role plot FAILED for %s: %s\n", cname, e$message))
    cat("         Typically means one group has too few pathways for 2D scatter.\n")
    try(safe_dev_off(), silent = TRUE)
    plot_section_failed <<- plot_section_failed + 1L
    role_plot_failed <<- role_plot_failed + 1L
  })
}



cat(sprintf("\n[summary] Signaling role plots: %d/%d succeeded, %d failed\n",
            role_plot_total - role_plot_failed, role_plot_total,
            role_plot_failed))

# ============================================================================
# 9. Export summary tables (controlled by save_qc_tables)
# ============================================================================
cat("\n[step] Exporting summary tables...\n")

if (isTRUE(cfg$output$save_qc_tables)) {
  # 9a. Per-group pathway list
  for (grp in analysis_groups) {
    pw <- cc_list[[grp]]@netP$pathways %||% character(0)
    if (length(pw) == 0) {
      dt <- data.table(group = character(0), pathway = character(0), rank = integer(0))
    } else {
      dt <- data.table(group = grp, pathway = pw, rank = seq_along(pw))
    }
    fwrite(dt, file.path(out$reports, sprintf("pathways_%s.csv", grp)))
  }

  # 9b. Pathway presence matrix
  all_pw <- sort(unique(unlist(lapply(cc_list, function(x) x@netP$pathways %||% character(0)))))
  pw_mat <- data.table(pathway = all_pw)

  for (grp in analysis_groups) {
    grp_pw <- cc_list[[grp]]@netP$pathways %||% character(0)
    pw_mat[[grp]] <- as.integer(all_pw %in% grp_pw)
  }

  if (nrow(pw_mat) > 0) {
    pw_mat[, n_groups := rowSums(.SD), .SDcols = analysis_groups]
    pw_mat <- pw_mat[order(-n_groups, pathway)]
  } else {
    pw_mat[, n_groups := integer(0)]
  }

  fwrite(pw_mat, file.path(out$reports, "pathway_presence_matrix.csv"))

  shared_all <- if (nrow(pw_mat) == 0) 0 else sum(pw_mat$n_groups == length(analysis_groups))
  unique_one <- if (nrow(pw_mat) == 0) 0 else sum(pw_mat$n_groups == 1)
  cat(sprintf("[report] %d total pathways: %d shared by all groups, %d unique to one group\n",
              nrow(pw_mat), shared_all, unique_one))

  # 9c. Per-group sender-receiver pair summary
  for (grp in analysis_groups) {
    net <- cc_list[[grp]]@net
    count_mat <- net$count
    weight_mat <- net$weight

    if (is.null(count_mat) || length(count_mat) == 0) {
      pairs <- data.table(sender = character(0), receiver = character(0),
                          n_interactions = numeric(0), interaction_weight = numeric(0))
    } else {
      ct_names <- rownames(count_mat)
      if (is.null(ct_names)) {
        ct_names <- as.character(seq_len(nrow(count_mat)))
      }

      pairs <- as.data.table(expand.grid(sender = ct_names, receiver = ct_names, stringsAsFactors = FALSE))
      pairs[, n_interactions := as.vector(count_mat)]
      pairs[, interaction_weight := as.vector(weight_mat)]
      pairs <- pairs[n_interactions > 0]
      setorder(pairs, -n_interactions, -interaction_weight)
    }

    fwrite(pairs, file.path(out$reports, sprintf("celltype_pairs_%s.csv", grp)))
  }
  cat("[report] celltype pair tables exported\n")

  # 9d. Pairwise pathway diff tables
  for (comp in comparisons) {
    cname <- comp$name
    num <- comp$numerator
    den <- comp$denominator

    pw_num <- cc_list[[num]]@netP$pathways %||% character(0)
    pw_den <- cc_list[[den]]@netP$pathways %||% character(0)

    report <- data.table(
      contrast = cname,
      contrast_label = label_contrast(cname, num, den),
      numerator_group = num,
      numerator_label = label_group(num),
      denominator_group = den,
      denominator_label = label_group(den),
      category = c(
        rep("num_only", length(setdiff(pw_num, pw_den))),
        rep("den_only", length(setdiff(pw_den, pw_num))),
        rep("shared", length(intersect(pw_num, pw_den)))
      ),
      pathway = c(
        setdiff(pw_num, pw_den),
        setdiff(pw_den, pw_num),
        intersect(pw_num, pw_den)
      )
    )

    fwrite(report, file.path(out$reports, sprintf("pathway_diff_%s.csv", cname)))
    cat(sprintf("  %s: %d shared, %d %s-only, %d %s-only\n",
                cname,
                sum(report$category == "shared"),
                sum(report$category == "num_only"), num,
                sum(report$category == "den_only"), den))
  }
} else {
  cat("[report] save_qc_tables=false -> summary CSV exports skipped\n")
}

# ============================================================================
# 10. Provenance
# ============================================================================
if (isTRUE(cfg$output$save_config_copy)) {
  invisible(file.copy(config_path, file.path(out$qc, "cellchat_params_snapshot_03.yaml"), overwrite = TRUE))
}

if (isTRUE(cfg$output$save_session_info)) {
  writeLines(capture.output(sessionInfo()), file.path(out$qc, "03_compare_session_info.txt"))
}

if (isTRUE(cfg$output$save_cellchat_version)) {
  writeLines(
    c(
      sprintf("CellChat: %s", packageVersion("CellChat")),
      sprintf("R:        %s", R.version.string),
      sprintf("Dataset:  %s", dataset_id),
      sprintf("Step:     03_compare_groups"),
      sprintf("Date:     %s", Sys.time())
    ),
    file.path(out$qc, "versions_03.txt")
  )
}

# ============================================================================
# 11. Summary
# ============================================================================
elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
comparison_names <- vapply(comparisons, `[[`, "", "name")

if (!isTRUE(cfg$output$save_qc_tables)) {
  cat("[done] save_qc_tables=false: reports are intentionally skipped\n")
}

# ★ FIX 5: Completion quality check — report any failures before done banner
if (plot_section_failed > 0L) {
  cat(sprintf(
    "\n[WARN] %d/%d plot section(s) failed during execution — review [WARN]/[SKIP] markers above.\n",
    plot_section_failed, plot_section_total
  ))
  cat("       Outputs are PARTIAL. Check logs for details.\n")
}

cat("\n")
rule <- strrep("=", 64)
cat(sprintf("+%s+\n", rule))
cat("| CellChat Step 03 Complete\n")
cat(sprintf("+%s+\n", rule))
cat(sprintf("| Dataset:      %s\n", dataset_id))
cat(sprintf("| Input mode:   %s\n", source_mode))
cat(sprintf("| Groups:       %s\n", paste(analysis_groups, collapse = ", ")))
cat(sprintf("| Comparisons:  %s\n", paste(comparison_names, collapse = ", ")))
cat(sprintf("| Plot issues:  %d/%d sections had errors\n", plot_section_failed, plot_section_total))
cat(sprintf("| Plots:        %s\n", out$plots))
cat(sprintf("| Reports:      %s\n", out$reports))
cat(sprintf("| Time:         %.1f min\n", elapsed))

next_script <- file.path(module_root, "scripts", "04_deep_dive.R")
if (file.exists(next_script)) {
  cat("| Next:         Rscript scripts/04_deep_dive.R\n")
} else {
  cat("| Next:         04_deep_dive.R (not yet created)\n")
}

cat(sprintf("+%s+\n", rule))
