# ============================================================================
# 04_deep_dive.R -- Per-pathway detailed CellChat visualization
#
# Input:    results/scrna/{dataset_id}/09_cellchat/objects/
# Output:   results/scrna/{dataset_id}/09_cellchat/deep_dive/{plots,reports}
# Usage:    cd modules/cellchat && Rscript scripts/04_deep_dive.R
# Override: DATASET_ID=gse253718 Rscript scripts/04_deep_dive.R
# ============================================================================

cat(
  "\n+================================================================+\n",
  "|             CellChat Step 04: Deep-Dive Analysis             |\n",
  "+================================================================+\n",
  sep = ""
)

suppressPackageStartupMessages({
  library(CellChat)
  library(yaml)
  library(here)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(dplyr)
  library(tidyr)
})

required_pkgs <- c("CellChat", "yaml", "here", "ggplot2", "patchwork", "data.table", "dplyr", "tidyr")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(sprintf(
    "[PREFLIGHT] Missing packages: %s\n  Fix: cd modules/cellchat && Rscript --vanilla -e 'renv::hydrate(); renv::snapshot()'",
    paste(missing_pkgs, collapse = ", ")
  ))
}

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

validate_group_membership <- function(analysis_groups, comparisons) {
  bad <- character(0)
  for (comp in comparisons) {
    if (!(comp$numerator %in% analysis_groups) || !(comp$denominator %in% analysis_groups)) {
      bad <- c(bad, sprintf("%s(%s vs %s)", comp$name, comp$numerator, comp$denominator))
    }
  }
  if (length(bad) > 0) {
    stop(sprintf(
      "[CONTRACT] Comparisons reference groups outside analysis scope: %s\n  analysis_groups: %s",
      paste(bad, collapse = ", "),
      paste(analysis_groups, collapse = ", ")
    ))
  }
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
      return(list(objects = cc_list[analysis_groups], source_mode = "cellchat_list_inferred"))
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
    return(list(objects = cc_list, source_mode = "per_group_inferred"))
  }

  per_group_legacy <- file.path(obj_dir, sprintf("cellchat_%s.rds", analysis_groups))
  if (all(file.exists(per_group_legacy))) {
    cc_list <- lapply(per_group_legacy, readRDS)
    names(cc_list) <- analysis_groups
    return(list(objects = cc_list, source_mode = "per_group_legacy"))
  }

  missing_inferred <- analysis_groups[!file.exists(per_group_inferred)]
  missing_legacy <- analysis_groups[!file.exists(per_group_legacy)]

  stop(sprintf(
    paste0(
      "[CONTRACT] Cannot resolve Step 04 input objects for all analysis groups.\n",
      "  Tried 1) cellchat_list_inferred.rds\n",
      "  Tried 2) cellchat_{group}_inferred.rds (missing: %s)\n",
      "  Tried 3) cellchat_{group}.rds (missing: %s)"
    ),
    if (length(missing_inferred) == 0) "none" else paste(missing_inferred, collapse = ", "),
    if (length(missing_legacy) == 0) "none" else paste(missing_legacy, collapse = ", ")
  ))
}

resolve_celltype_indices <- function(pattern, all_ct) {
  exact <- which(all_ct == pattern)
  if (length(exact) > 0) return(exact)
  partial <- grep(pattern, all_ct, ignore.case = TRUE)
  if (length(partial) > 0) return(partial)
  integer(0)
}

adaptive_net_size <- function(base_size, n_celltypes) {
  max(base_size, base_size + (n_celltypes - 10) * 0.2)
}

t_start <- Sys.time()
plot_section_total <- 0L
plot_section_failed <- 0L

# ============================================================================
# 1. Bootstrap config and output paths
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
  deep = file.path(out_base, "deep_dive"),
  deep_plots = file.path(out_base, "deep_dive", "plots"),
  deep_reports = file.path(out_base, "deep_dive", "reports")
)

for (d in out[c("objects", "qc", "deep", "deep_plots", "deep_reports")]) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat(sprintf("[config] Loaded: %s\n", config_path))
cat(sprintf("[dataset] id: %s\n", dataset_id))
cat(sprintf("[paths] module_root:  %s\n", module_root))
cat(sprintf("[paths] project_root: %s\n", project_root))
cat(sprintf("[paths] out_base:     %s\n", out$base))
cat(sprintf("[env] CellChat %s\n", as.character(packageVersion("CellChat"))))

# ============================================================================
# 2. Load Step01 design contract
# ============================================================================
design_path <- file.path(out$objects, "design.rds")
if (!file.exists(design_path)) {
  stop("[CONTRACT] design.rds not found. Run 01_prepare_cellchat.R first.\n  Expected: ", design_path)
}
design <- readRDS(design_path)

if (is.null(design$all_groups) || length(design$all_groups) == 0) {
  stop("[CONTRACT] design.rds missing non-empty all_groups")
}
if (is.null(design$comparisons) || length(design$comparisons) == 0) {
  stop("[CONTRACT] design.rds missing non-empty comparisons")
}

analysis_groups <- unique(as.character(design$all_groups))
reference_group <- design$reference_group %||% NULL
if (!is.null(reference_group) && reference_group %in% analysis_groups) {
  analysis_groups <- c(reference_group, setdiff(analysis_groups, reference_group))
}
comparisons <- normalize_comparisons(design$comparisons)
validate_group_membership(analysis_groups, comparisons)

cat(sprintf("[load] Design groups (%d): %s\n", length(analysis_groups), paste(analysis_groups, collapse = ", ")))
cat(sprintf("[load] Comparisons (%d): %s\n",
            length(comparisons), paste(vapply(comparisons, `[[`, "", "name"), collapse = ", ")))

# ============================================================================
# 3. Deep-dive config
# ============================================================================
comparison_cfg <- cfg$comparison %||% list()
dd_cfg <- cfg$deep_dive %||% list()

dd_pathways <- dd_cfg$pathways
if (is.null(dd_pathways) || length(dd_pathways) == 0) {
  dd_pathways <- comparison_cfg$focus_pathways
}
dd_pathways <- unique(as.character(dd_pathways %||% character(0)))
if (length(dd_pathways) == 0) {
  stop("[CONTRACT] No pathways configured for deep-dive. Set deep_dive.pathways or comparison.focus_pathways")
}

ct_pairs <- dd_cfg$cell_type_pairs %||% list()
if (!is.list(ct_pairs)) {
  stop("[CONTRACT] deep_dive.cell_type_pairs must be a list")
}
invalid_pairs <- which(!vapply(ct_pairs, function(x) is.list(x) && !is.null(x$source) && !is.null(x$target), logical(1)))
if (length(invalid_pairs) > 0) {
  stop(sprintf("[CONTRACT] Invalid deep_dive.cell_type_pairs entries at indices: %s",
               paste(invalid_pairs, collapse = ", ")))
}

dd_layouts <- as.character(dd_cfg$layouts %||% c("chord", "circle"))
valid_layouts <- c("chord", "circle")
bad_layouts <- setdiff(dd_layouts, valid_layouts)
if (length(bad_layouts) > 0) {
  stop(sprintf("[CONTRACT] Unsupported deep_dive.layouts values: %s", paste(bad_layouts, collapse = ", ")))
}

chord_size <- as.numeric(dd_cfg$chord_size %||% 10)
circle_size <- as.numeric(dd_cfg$circle_size %||% 10)
contribution_top_n <- as.integer(dd_cfg$contribution_top_n %||% 20)
heatmap_font_size <- as.numeric(dd_cfg$heatmap_font_size %||% 9)
heatmap_color <- as.character(dd_cfg$heatmap_color %||% "Reds")
pair_bubble_font <- as.numeric(dd_cfg$pair_bubble_font_size %||% 10)
pair_bubble_thresh <- dd_cfg$pair_bubble_thresh
if (!is.null(pair_bubble_thresh) && (!is.numeric(pair_bubble_thresh) || length(pair_bubble_thresh) != 1 || is.na(pair_bubble_thresh))) {
  stop("[CONTRACT] deep_dive.pair_bubble_thresh must be numeric scalar or null")
}

if (!is.finite(chord_size) || chord_size <= 0) stop("[CONTRACT] deep_dive.chord_size must be > 0")
if (!is.finite(circle_size) || circle_size <= 0) stop("[CONTRACT] deep_dive.circle_size must be > 0")
if (!is.finite(heatmap_font_size) || heatmap_font_size <= 0) stop("[CONTRACT] deep_dive.heatmap_font_size must be > 0")
if (!is.finite(pair_bubble_font) || pair_bubble_font <= 0) stop("[CONTRACT] deep_dive.pair_bubble_font_size must be > 0")
if (is.na(contribution_top_n) || contribution_top_n < 1) stop("[CONTRACT] deep_dive.contribution_top_n must be >= 1")

cat(sprintf("[config] deep_dive.pathways: %d\n", length(dd_pathways)))
cat(sprintf("[config] deep_dive.cell_type_pairs: %d\n", length(ct_pairs)))
cat(sprintf("[config] deep_dive.layouts: %s\n", paste(dd_layouts, collapse = ", ")))

# ============================================================================
# 4. Load objects from Step02/03 artifacts
# ============================================================================
loaded <- load_cellchat_objects(out$objects, analysis_groups)
cc_list <- loaded$objects
source_mode <- loaded$source_mode
cat(sprintf("[load] Input mode: %s\n", source_mode))

for (grp in analysis_groups) {
  cc <- cc_list[[grp]]
  if (is.null(cc@idents) || length(cc@idents) == 0) {
    stop(sprintf("[CONTRACT] %s object missing @idents.", grp))
  }
  if (is.null(cc@net) || length(cc@net) == 0 || is.null(cc@netP) || length(cc@netP) == 0) {
    stop(sprintf("[CONTRACT] %s object missing inferred slots (@net/@netP).", grp))
  }
  cat(sprintf("[load] %s: %d cells, %d pathways\n",
              grp, length(cc@idents), length(cc@netP$pathways %||% character(0))))
}

merged_path <- file.path(out$objects, "cellchat_merged.rds")
if (!file.exists(merged_path)) {
  stop("[CONTRACT] Missing merged object: ", merged_path)
}
cc_merged <- readRDS(merged_path)

if (is.null(cc_merged@meta) || !"datasets" %in% colnames(cc_merged@meta)) {
  stop("[CONTRACT] merged object missing meta$datasets")
}
merged_groups <- unique(as.character(cc_merged@meta$datasets))
missing_merged <- setdiff(analysis_groups, merged_groups)
if (length(missing_merged) > 0) {
  stop(sprintf("[CONTRACT] merged object missing analysis groups: %s",
               paste(missing_merged, collapse = ", ")))
}
cat(sprintf("[load] Merged datasets: %d (%s)\n", length(merged_groups), paste(merged_groups, collapse = ", ")))

all_celltypes <- levels(cc_list[[analysis_groups[1]]]@idents)
n_ct <- length(all_celltypes)
cat(sprintf("[load] Cell types: %d\n", n_ct))

# ============================================================================
# 5. Per-pathway deep-dive
# ============================================================================
cat("\n[step] Per-pathway deep-dive analysis...\n")
cat(sprintf("[step] Pathways: %d | Comparisons: %d | Layouts: %s\n",
            length(dd_pathways), length(comparisons), paste(dd_layouts, collapse = "+")))

pw_summary_rows <- list()

for (pw in dd_pathways) {
  cat(sprintf("\n  === Pathway: %s ===\n", pw))

  pw_in_group <- sapply(analysis_groups, function(g) {
    pw %in% (cc_list[[g]]@netP$pathways %||% character(0))
  })
  groups_with_pw <- analysis_groups[pw_in_group]
  groups_without_pw <- analysis_groups[!pw_in_group]

  if (length(groups_without_pw) > 0) {
    cat(sprintf("  [note] Not detected in: %s\n", paste(groups_without_pw, collapse = ", ")))
  }
  if (length(groups_with_pw) == 0) {
    cat(sprintf("  [skip] %s not detected in any group\n", pw))
    next
  }

  pw_safe <- gsub("[^A-Za-z0-9_-]", "_", pw)
  n_panels <- length(groups_with_pw)
  grp_tag <- paste(gsub("[^A-Za-z0-9_-]", "", groups_with_pw), collapse = "_")
  if (nchar(grp_tag) == 0) grp_tag <- "none"

  if ("chord" %in% dd_layouts) {
    plot_section_total <- plot_section_total + 1L
    tryCatch({
      sz <- adaptive_net_size(chord_size, n_ct)
      chord_stem <- sprintf("01_chord_%s_%dgrp_%s", pw_safe, n_panels, grp_tag)

      pdf(file.path(out$deep_plots, paste0(chord_stem, ".pdf")),
          width = sz * n_panels, height = sz)
      par(mfrow = c(1, n_panels), mar = c(1, 1, 3, 1))
      for (grp in groups_with_pw) {
        netVisual_aggregate(
          cc_list[[grp]],
          signaling = pw,
          layout = "chord",
          title.name = sprintf("%s - %s", pw, label_group(grp)),
          show.legend = (grp == groups_with_pw[length(groups_with_pw)])
        )
      }
      safe_dev_off()

      png(file.path(out$deep_plots, paste0(chord_stem, ".png")),
          width = sz * n_panels * 100, height = sz * 100, res = 150)
      par(mfrow = c(1, n_panels), mar = c(1, 1, 3, 1))
      for (grp in groups_with_pw) {
        netVisual_aggregate(
          cc_list[[grp]],
          signaling = pw,
          layout = "chord",
          title.name = sprintf("%s - %s", pw, label_group(grp)),
          show.legend = (grp == groups_with_pw[length(groups_with_pw)])
        )
      }
      safe_dev_off()

      cat(sprintf("  [plot] %s (%d panels)\n", chord_stem, n_panels))
    }, error = function(e) {
      cat(sprintf("  [WARN] Chord failed for %s: %s\n", pw, e$message))
      try(safe_dev_off(), silent = TRUE)
      plot_section_failed <<- plot_section_failed + 1L
    })
  }

  if ("circle" %in% dd_layouts) {
    plot_section_total <- plot_section_total + 1L
    tryCatch({
      sz <- adaptive_net_size(circle_size, n_ct)
      circle_stem <- sprintf("02_circle_%s_%dgrp_%s", pw_safe, n_panels, grp_tag)

      pdf(file.path(out$deep_plots, paste0(circle_stem, ".pdf")),
          width = sz * n_panels, height = sz)
      par(mfrow = c(1, n_panels), mar = c(1, 1, 3, 1))
      for (grp in groups_with_pw) {
        netVisual_aggregate(
          cc_list[[grp]],
          signaling = pw,
          layout = "circle"
        )
        title(main = sprintf("%s - %s", pw, label_group(grp)), cex.main = 0.9, line = 0.5)
      }
      safe_dev_off()

      png(file.path(out$deep_plots, paste0(circle_stem, ".png")),
          width = sz * n_panels * 100, height = sz * 100, res = 150)
      par(mfrow = c(1, n_panels), mar = c(1, 1, 3, 1))
      for (grp in groups_with_pw) {
        netVisual_aggregate(
          cc_list[[grp]],
          signaling = pw,
          layout = "circle"
        )
        title(main = sprintf("%s - %s", pw, label_group(grp)), cex.main = 0.9, line = 0.5)
      }
      safe_dev_off()

      cat(sprintf("  [plot] %s (%d panels)\n", circle_stem, n_panels))
    }, error = function(e) {
      cat(sprintf("  [WARN] Circle failed for %s: %s\n", pw, e$message))
      try(safe_dev_off(), silent = TRUE)
      plot_section_failed <<- plot_section_failed + 1L
    })
  }

  plot_section_total <- plot_section_total + 1L
  tryCatch({
    lr_plots <- list()
    lr_tables <- list()
    for (grp in groups_with_pw) {
      p_lr <- netAnalysis_contribution(
        cc_list[[grp]],
        signaling = pw,
        title = sprintf("%s - %s", pw, grp)
      )
      lr_plots[[grp]] <- p_lr

      tryCatch({
        lr_tbl <- p_lr$data
        if (!is.null(lr_tbl)) {
          lr_tbl$group <- grp
          lr_tbl$group_label <- label_group(grp)
          lr_tbl$pathway <- pw
          lr_tables[[grp]] <- lr_tbl
        }
      }, error = function(e) NULL)
    }

    if (length(lr_plots) > 0) {
      n_lr <- length(lr_plots)
      combined_lr <- wrap_plots(lr_plots, ncol = n_lr)
      w_lr <- max(7, 6 * n_lr)
      h_lr <- max(6, min(12, contribution_top_n * 0.3 + 2))
      lr_stem <- sprintf("03_lr_contribution_%s_%dgrp_%s", pw_safe, n_lr, grp_tag)

      pdf(file.path(out$deep_plots, paste0(lr_stem, ".pdf")),
          width = w_lr, height = h_lr)
      print(combined_lr)
      safe_dev_off()

      png(file.path(out$deep_plots, paste0(lr_stem, ".png")),
          width = w_lr * 100, height = h_lr * 100, res = 150)
      print(combined_lr)
      safe_dev_off()

      if (isTRUE(cfg$output$save_qc_tables) && length(lr_tables) > 0) {
        lr_combined <- bind_rows(lr_tables)
        fwrite(lr_combined, file.path(out$deep_reports, paste0(lr_stem, ".csv")))
      }

      cat(sprintf("  [plot] %s (%d groups)\n", lr_stem, n_lr))
    }
  }, error = function(e) {
    cat(sprintf("  [WARN] LR contribution failed for %s: %s\n", pw, e$message))
    try(safe_dev_off(), silent = TRUE)
    plot_section_failed <<- plot_section_failed + 1L
  })

  plot_section_total <- plot_section_total + 1L
  tryCatch({
    sz_hm <- max(8, n_ct * 0.4 + 2)
    hm_stem <- sprintf("04_heatmap_%s_%dgrp_%s", pw_safe, n_panels, grp_tag)

    pdf(file.path(out$deep_plots, paste0(hm_stem, ".pdf")),
        width = sz_hm * n_panels + 2, height = sz_hm)
    for (grp in groups_with_pw) {
      tryCatch({
        netVisual_heatmap(
          cc_list[[grp]],
          signaling = pw,
          color.heatmap = heatmap_color,
          font.size = heatmap_font_size,
          title.name = sprintf("%s - %s", pw, label_group(grp))
        )
      }, error = function(e) {
        cat(sprintf("    [note] Heatmap skip %s/%s: %s\n", pw, grp, e$message))
      })
    }
    safe_dev_off()

    png(file.path(out$deep_plots, paste0(hm_stem, ".png")),
        width = (sz_hm * n_panels + 2) * 100, height = sz_hm * 100, res = 150)
    for (grp in groups_with_pw) {
      tryCatch({
        netVisual_heatmap(
          cc_list[[grp]],
          signaling = pw,
          color.heatmap = heatmap_color,
          font.size = heatmap_font_size,
          title.name = sprintf("%s - %s", pw, label_group(grp))
        )
      }, error = function(e) NULL)
    }
    safe_dev_off()

    cat(sprintf("  [plot] %s (%d panels)\n", hm_stem, n_panels))
  }, error = function(e) {
    cat(sprintf("  [WARN] Heatmap failed for %s: %s\n", pw, e$message))
    try(safe_dev_off(), silent = TRUE)
    plot_section_failed <<- plot_section_failed + 1L
  })

  plot_section_total <- plot_section_total + 1L
  tryCatch({
    sz_gc <- adaptive_net_size(chord_size, n_ct)
    chord_gene_stem <- sprintf("05_chord_gene_%s_%dgrp_%s", pw_safe, n_panels, grp_tag)

    pdf(file.path(out$deep_plots, paste0(chord_gene_stem, ".pdf")),
        width = sz_gc * n_panels, height = sz_gc)
    par(mfrow = c(1, n_panels), mar = c(1, 1, 3, 1))
    for (grp in groups_with_pw) {
      tryCatch({
        netVisual_chord_gene(
          cc_list[[grp]],
          signaling = pw,
          title.name = sprintf("%s genes - %s", pw, label_group(grp)),
          legend.pos.x = 5
        )
      }, error = function(e) {
        plot.new()
        title(sprintf("%s/%s: %s", pw, label_group(grp), e$message), cex.main = 0.8)
      })
    }
    safe_dev_off()

    png(file.path(out$deep_plots, paste0(chord_gene_stem, ".png")),
        width = sz_gc * n_panels * 100, height = sz_gc * 100, res = 150)
    par(mfrow = c(1, n_panels), mar = c(1, 1, 3, 1))
    for (grp in groups_with_pw) {
      tryCatch({
        netVisual_chord_gene(
          cc_list[[grp]],
          signaling = pw,
          title.name = sprintf("%s genes - %s", pw, label_group(grp)),
          legend.pos.x = 5
        )
      }, error = function(e) {
        plot.new()
        title(sprintf("%s/%s: skipped", pw, label_group(grp)), cex.main = 0.8)
      })
    }
    safe_dev_off()

    cat(sprintf("  [plot] %s (%d panels)\n", chord_gene_stem, n_panels))
  }, error = function(e) {
    cat(sprintf("  [WARN] Chord gene failed for %s: %s\n", pw, e$message))
    try(safe_dev_off(), silent = TRUE)
    plot_section_failed <<- plot_section_failed + 1L
  })

  for (grp in analysis_groups) {
    pres <- pw %in% (cc_list[[grp]]@netP$pathways %||% character(0))
    prob_val <- NA_real_
    if (pres) {
      tryCatch({
        prob_val <- sum(cc_list[[grp]]@netP$prob[pw, , ], na.rm = TRUE)
      }, error = function(e) NULL)
    }
    pw_summary_rows[[length(pw_summary_rows) + 1]] <- data.frame(
      pathway = pw,
      group = grp,
      detected = pres,
      total_prob = prob_val,
      stringsAsFactors = FALSE
    )
  }
}

# ============================================================================
# 6. Cell type pair focused analysis
# ============================================================================
cat("\n[step] Cell type pair focused analysis...\n")

if (length(ct_pairs) == 0) {
  cat("[skip] No cell_type_pairs configured - skipping pair-focused section\n")
} else {
  for (pi in seq_along(ct_pairs)) {
    pair <- ct_pairs[[pi]]
    src_name <- as.character(pair$source)
    tgt_name <- as.character(pair$target)

    src_idx <- resolve_celltype_indices(src_name, all_celltypes)
    tgt_idx <- resolve_celltype_indices(tgt_name, all_celltypes)

    if (length(src_idx) == 0 || length(tgt_idx) == 0) {
      cat(sprintf("[skip] Pair %d: '%s' -> '%s' not found in celltypes\n", pi, src_name, tgt_name))
      next
    }

    src_resolved <- all_celltypes[src_idx]
    tgt_resolved <- all_celltypes[tgt_idx]
    pair_label <- sprintf("%s_to_%s", gsub("[^A-Za-z0-9]", "", src_name), gsub("[^A-Za-z0-9]", "", tgt_name))

    cat(sprintf("[pair] %d: %s -> %s\n", pi,
                paste(src_resolved, collapse = "+"),
                paste(tgt_resolved, collapse = "+")))

    for (comp in comparisons) {
      cname <- comp$name
      num <- comp$numerator
      den <- comp$denominator

      idx_num <- match(num, analysis_groups)
      idx_den <- match(den, analysis_groups)
      if (is.na(idx_num) || is.na(idx_den)) next

      plot_section_total <- plot_section_total + 1L
      tryCatch({
        bubble_args <- list(
          object = cc_merged,
          comparison = c(idx_den, idx_num),
          signaling = dd_pathways,
          sources.use = src_idx,
          targets.use = tgt_idx,
          angle.x = 45,
          remove.isolate = TRUE,
          font.size = pair_bubble_font,
          title.name = sprintf("%s -> %s : %s",
                               paste(src_resolved, collapse = "/"),
                               paste(tgt_resolved, collapse = "/"),
                               label_contrast(cname, num, den))
        )
        if (!is.null(pair_bubble_thresh)) bubble_args$thresh <- pair_bubble_thresh

        p_pair <- do.call(netVisual_bubble, bubble_args)

        n_y_est <- max(1, length(dd_pathways)) * 3
        h_pair <- max(7, n_y_est * 0.3 + 2)
        w_pair <- 14

        pdf(file.path(out$deep_plots, sprintf("06_pair_bubble_%s_%s.pdf", pair_label, cname)),
            width = w_pair, height = h_pair)
        print(p_pair)
        safe_dev_off()

        png(file.path(out$deep_plots, sprintf("06_pair_bubble_%s_%s.png", pair_label, cname)),
            width = w_pair * 100, height = h_pair * 100, res = 150)
        print(p_pair)
        safe_dev_off()

        cat(sprintf("  [plot] 06_pair_bubble_%s_%s\n", pair_label, cname))
      }, error = function(e) {
        cat(sprintf("  [WARN] Pair bubble failed %s/%s: %s\n", pair_label, cname, e$message))
        try(safe_dev_off(), silent = TRUE)
        plot_section_failed <<- plot_section_failed + 1L
      })
    }

    plot_section_total <- plot_section_total + 1L
    tryCatch({
      prob_rows <- list()
      for (grp in analysis_groups) {
        cc <- cc_list[[grp]]
        prob_mat <- cc@net$prob
        if (is.null(prob_mat)) next

        for (si in src_idx) {
          for (ti in tgt_idx) {
            total_p <- sum(prob_mat[si, ti, ], na.rm = TRUE)
            prob_rows[[length(prob_rows) + 1]] <- data.frame(
              source = all_celltypes[si],
              target = all_celltypes[ti],
              group = grp,
              group_label = label_group(grp),
              total_prob = total_p,
              stringsAsFactors = FALSE
            )
          }
        }
      }

      if (length(prob_rows) > 0) {
        prob_df <- bind_rows(prob_rows)

        p_bar <- ggplot(prob_df, aes(x = group_label, y = total_prob, fill = group_label)) +
          geom_col(width = 0.6) +
          facet_wrap(~ paste(source, "->", target), scales = "free_y") +
          labs(
            title = sprintf("Communication probability: %s -> %s",
                            paste(src_resolved, collapse = "/"),
                            paste(tgt_resolved, collapse = "/")),
            subtitle = "Control (C) vs Xan (K)",
            x = "Group",
            y = "Total probability (sum over all LR pairs)"
          ) +
          theme_minimal(base_size = 12) +
          theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 30, hjust = 1)
          )

        n_facets <- nrow(distinct(prob_df, source, target))
        w_bar <- max(8, n_facets * 4)
        h_bar <- 6

        pdf(file.path(out$deep_plots, sprintf("07_pair_prob_%s.pdf", pair_label)),
            width = w_bar, height = h_bar)
        print(p_bar)
        safe_dev_off()

        png(file.path(out$deep_plots, sprintf("07_pair_prob_%s.png", pair_label)),
            width = w_bar * 100, height = h_bar * 100, res = 150)
        print(p_bar)
        safe_dev_off()

        if (isTRUE(cfg$output$save_qc_tables)) {
          fwrite(prob_df, file.path(out$deep_reports, sprintf("07_pair_prob_%s.csv", pair_label)))
        }

        cat(sprintf("  [plot] 07_pair_prob_%s (%d facets)\n", pair_label, n_facets))
      }
    }, error = function(e) {
      cat(sprintf("  [WARN] Pair prob failed %s: %s\n", pair_label, e$message))
      try(safe_dev_off(), silent = TRUE)
      plot_section_failed <<- plot_section_failed + 1L
    })
  }
}

# ============================================================================
# 7. Reports / provenance
# ============================================================================
cat("\n[step] Writing summary reports...\n")

if (isTRUE(cfg$output$save_qc_tables) && length(pw_summary_rows) > 0) {
  pw_summary <- bind_rows(pw_summary_rows)
  pw_summary$group_label <- label_group(pw_summary$group)
  fwrite(pw_summary, file.path(out$deep_reports, "pathway_group_summary.csv"))
  cat(sprintf("[report] pathway_group_summary.csv: %d rows\n", nrow(pw_summary)))

  pw_wide <- pw_summary %>%
    select(pathway, group, total_prob) %>%
    tidyr::pivot_wider(names_from = group, values_from = total_prob)
  fwrite(pw_wide, file.path(out$deep_reports, "pathway_group_probability_wide.csv"))
  cat(sprintf("[report] pathway_group_probability_wide.csv: %d pathways\n", nrow(pw_wide)))
} else if (!isTRUE(cfg$output$save_qc_tables)) {
  cat("[report] save_qc_tables=false -> deep-dive CSV exports skipped\n")
}

if (isTRUE(cfg$output$save_config_copy)) {
  invisible(file.copy(config_path, file.path(out$qc, "cellchat_params_snapshot_04.yaml"), overwrite = TRUE))
}

if (isTRUE(cfg$output$save_session_info)) {
  writeLines(capture.output(sessionInfo()), file.path(out$deep_reports, "sessionInfo_step04.txt"))
}

if (isTRUE(cfg$output$save_cellchat_version)) {
  writeLines(
    c(
      sprintf("CellChat: %s", packageVersion("CellChat")),
      sprintf("R:        %s", R.version.string),
      sprintf("Dataset:  %s", dataset_id),
      sprintf("Step:     04_deep_dive"),
      sprintf("Date:     %s", Sys.time())
    ),
    file.path(out$qc, "versions_04.txt")
  )
}

all_files <- list.files(c(out$deep_plots, out$deep_reports), full.names = TRUE, recursive = TRUE)
if (length(all_files) > 0) {
  all_info <- file.info(all_files)
  run_files <- all_files[!is.na(all_info$mtime) & all_info$mtime >= (t_start - 1)]
} else {
  run_files <- character(0)
}
cat(sprintf("[report] Output files (this run): %d | total present: %d\n",
            length(run_files), length(all_files)))

# ============================================================================
# 8. Final summary
# ============================================================================
elapsed <- round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1)

cat("\n")
rule <- strrep("=", 64)
cat(sprintf("+%s+\n", rule))
cat("| CellChat Step 04 Complete\n")
cat(sprintf("+%s+\n", rule))
cat(sprintf("| Dataset:      %s\n", dataset_id))
cat(sprintf("| Input mode:   %s\n", source_mode))
cat(sprintf("| Pathways:     %d\n", length(dd_pathways)))
cat(sprintf("| CT pairs:     %d\n", length(ct_pairs)))
cat(sprintf("| Comparisons:  %d\n", length(comparisons)))
cat(sprintf("| Plot issues:  %d/%d sections had errors\n", plot_section_failed, plot_section_total))
cat(sprintf("| Plots:        %s\n", out$deep_plots))
cat(sprintf("| Reports:      %s\n", out$deep_reports))
cat(sprintf("| Time:         %.1f min\n", elapsed))
cat(sprintf("+%s+\n", rule))
