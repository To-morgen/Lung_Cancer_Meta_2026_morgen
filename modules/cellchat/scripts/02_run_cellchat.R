# ============================================================================
# 02_run_cellchat.R — Run CellChat inference per group
#
# Module:  modules/cellchat/  (stateless)
# Input:   objects/cellchat_{group}.rds  +  objects/design.rds
# Output:  objects/cellchat_{group}_inferred.rds  +  qc/  +  reports/
#
# Core pipeline per group:
#   subsetData → identifyOverExpressedGenes → identifyOverExpressedInteractions
#   → computeCommunProb → filterCommunication
#   → computeCommunProbPathway → aggregateNet
#
# Usage:   cd modules/cellchat && Rscript scripts/02_run_cellchat.R
# ============================================================================

cat("
╔══════════════════════════════════════════════════════════════╗
║       CellChat Step 02: Run Inference Per Group             ║
╚══════════════════════════════════════════════════════════════╝
\n")

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(CellChat)
  library(data.table)
})

# ============================================================================
# 0. Path resolution & config
# ============================================================================
module_root  <- here::here()
project_root <- normalizePath(file.path(module_root, "..", ".."))

config_path <- file.path(module_root, "configs", "cellchat_params.yaml")
if (!file.exists(config_path)) stop("[CONTRACT] Config not found: ", config_path)
cfg <- yaml::read_yaml(config_path)
cat(sprintf("[config] Loaded: %s\n", config_path))

# Dataset ID
ds_id <- Sys.getenv("DATASET_ID", "")
if (ds_id == "" && !is.null(cfg$dataset$id)) ds_id <- cfg$dataset$id
if (ds_id == "") stop("[CONTRACT] No dataset_id")
cat(sprintf("[dataset] id: %s\n", ds_id))

# Output base
base_resolved <- gsub("\\{dataset_id\\}", ds_id, cfg$output$base_dir)
out_base <- file.path(project_root, base_resolved)

out <- list(
  base    = out_base,
  objects = file.path(out_base, "objects"),
  qc      = file.path(out_base, "qc"),
  reports = file.path(out_base, "reports")
)
for (d in out) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 1. Load design + CellChat objects from Step 01
# ============================================================================
design_path <- file.path(out$objects, "design.rds")
if (!file.exists(design_path)) {
  stop("[CONTRACT] design.rds not found. Run 01_prepare_cellchat.R first.\n  Expected: ",
       design_path)
}
design <- readRDS(design_path)
cat(sprintf("[load] Design: %d groups, %d comparisons\n",
            length(design$all_groups), length(design$comparisons)))

list_path <- file.path(out$objects, "cellchat_list.rds")
if (!file.exists(list_path)) {
  stop("[CONTRACT] cellchat_list.rds not found. Run 01_prepare_cellchat.R first.")
}
cellchat_list <- readRDS(list_path)
cat(sprintf("[load] cellchat_list groups: %s\n", paste(names(cellchat_list), collapse = ", ")))

# ============================================================================
# 1b. Input contract validation
# ============================================================================
# Group coverage
missing_in_list <- setdiff(design$all_groups, names(cellchat_list))
extra_in_list   <- setdiff(names(cellchat_list), design$all_groups)
if (length(missing_in_list) > 0) {
  stop(sprintf("[CONTRACT] design$all_groups has groups not in cellchat_list: %s",
               paste(missing_in_list, collapse = ", ")))
}
if (length(extra_in_list) > 0) {
  cat(sprintf("[WARN] cellchat_list has extra groups not in design: %s (will be ignored)\n",
              paste(extra_in_list, collapse = ", ")))
}

analysis_groups <- design$all_groups
cat(sprintf("[validate] Analysis scope (design$all_groups): %s\n",
            paste(analysis_groups, collapse = ", ")))

# Per-object structural check
for (grp in analysis_groups) {
  cc <- cellchat_list[[grp]]
  if (is.null(cc@DB) || length(cc@DB) == 0) {
    stop(sprintf("[CONTRACT] cellchat_list[['%s']]@DB is empty. Step 01 did not assign CellChatDB.", grp))
  }
  if (is.null(cc@idents) || length(cc@idents) == 0) {
    stop(sprintf("[CONTRACT] cellchat_list[['%s']]@idents is empty. Step 01 object is malformed.", grp))
  }
  cat(sprintf("[load] %s: %d cells, %d celltypes\n",
              grp, length(cc@idents), length(levels(cc@idents))))
}
cat("[validate] All CellChat objects pass structural contract\n")

# ============================================================================
# 2. Inference parameters
# ============================================================================
inf <- cfg$inference
cat(sprintf("[params] type=%s, nboot=%d, thresh=%.3f, seed=%d\n",
            inf$type, inf$nboot, inf$thresh, inf$seed))

pop_size <- isTRUE(inf$population.size)
cat(sprintf("[params] population.size=%s\n", pop_size))

# population.size_sensitivity: RESERVED — not consumed by current CellChat v2 API
if (!is.null(inf$population.size_sensitivity)) {
  cat(sprintf("[params] population.size_sensitivity=%s (RESERVED — not consumed by current CellChat API)\n",
              inf$population.size_sensitivity))
}

# Validate trim
if (!is.null(inf$trim)) {
  if (!is.numeric(inf$trim) || inf$trim < 0 || inf$trim > 0.5) {
    stop(sprintf("[CONTRACT] inference.trim must be numeric in [0, 0.5], got: %s",
                 deparse(inf$trim)))
  }
  if (inf$type != "truncatedMean") {
    cat(sprintf("[WARN] trim=%.2f specified but type='%s' (trim only effective with 'truncatedMean')\n",
                inf$trim, inf$type))
  }
}

# ============================================================================
# 3. Run inference loop
# ============================================================================
inference_summary <- list()

for (grp in analysis_groups) {
  cat(sprintf("\n%s\n", strrep("━", 64)))
  cat(sprintf("  Running CellChat inference: %s\n", grp))
  cat(sprintf("%s\n", strrep("━", 64)))

  cc <- cellchat_list[[grp]]
  t_start <- Sys.time()

  # ── 3a. Preprocessing ──
  cat("[step] subsetData...\n")
  cc <- subsetData(cc)

  if (is.null(cc@data.signaling) || nrow(cc@data.signaling) == 0) {
    stop(sprintf(
      "[FATAL] %s: subsetData produced empty object@data.signaling. Check database.species='%s' and gene symbol format.",
      grp, cfg$database$species
    ))
  }

  cat(sprintf("[step] identifyOverExpressedGenes (%d signaling genes)...\n",
              nrow(cc@data.signaling)))
  cc <- identifyOverExpressedGenes(cc)

  cat("[step] identifyOverExpressedInteractions...\n")
  cc <- identifyOverExpressedInteractions(cc)

  # ── 3b. Communication probability ──
  cat(sprintf("[step] computeCommunProb (type=%s, nboot=%d, population.size=%s)...\n",
              inf$type, inf$nboot, pop_size))

  set.seed(inf$seed)
  cc <- computeCommunProb(
    cc,
    type            = inf$type,
    trim            = inf$trim,
    population.size = pop_size,
    nboot           = inf$nboot,
    seed.use        = inf$seed
  )

  # ── 3c. Filter low-confidence interactions ──
  cat(sprintf("[step] filterCommunication (min.cells=%d)...\n",
              cfg$filtering$min_cells_per_celltype))
  cc <- filterCommunication(cc, min.cells = cfg$filtering$min_cells_per_celltype)

  # ── 3d. Pathway-level aggregation ──
  cat(sprintf("[step] computeCommunProbPathway (thresh=%.3f)...\n", inf$thresh))
  cc <- computeCommunProbPathway(cc, thresh = inf$thresh)

  cat("[step] aggregateNet...\n")
  cc <- aggregateNet(cc)

  t_end <- Sys.time()
  elapsed <- round(as.numeric(difftime(t_end, t_start, units = "mins")), 1)

  # ── 3e. Quick summary ──
  n_LR <- nrow(subsetCommunication(cc))
  n_pathway <- length(cc@netP$pathways)

  cat(sprintf("[done] %s: %d L-R pairs, %d pathways, %.1f min\n",
              grp, n_LR, n_pathway, elapsed))

  inference_summary[[grp]] <- data.frame(
    group       = grp,
    n_cells     = length(cc@idents),
    n_celltypes = length(levels(cc@idents)),
    n_LR_pairs  = n_LR,
    n_pathways  = n_pathway,
    elapsed_min = elapsed,
    stringsAsFactors = FALSE
  )

  cellchat_list[[grp]] <- cc
}

# ============================================================================
# 4. Save inferred objects
# ============================================================================
cat("\n[save] Saving inferred CellChat objects...\n")

for (grp in analysis_groups) {
  fpath <- file.path(out$objects, sprintf("cellchat_%s_inferred.rds", grp))
  saveRDS(cellchat_list[[grp]], fpath)
  cat(sprintf("  %s\n", fpath))
}

# Save list with inferred versions
list_inferred_path <- file.path(out$objects, "cellchat_list_inferred.rds")
saveRDS(cellchat_list[analysis_groups], list_inferred_path)
cat(sprintf("  %s\n", list_inferred_path))

# ============================================================================
# 5. Export L-R communication tables
# ============================================================================
cat("\n[export] Exporting L-R communication tables...\n")

for (grp in analysis_groups) {
  cc <- cellchat_list[[grp]]

  # Full L-R table
  lr_df <- subsetCommunication(cc)
  if (nrow(lr_df) > 0) {
    fwrite(lr_df, file.path(out$reports, sprintf("LR_communication_%s.csv", grp)))
    cat(sprintf("  %s: %d L-R pairs\n", grp, nrow(lr_df)))
  } else {
    cat(sprintf("  %s: 0 L-R pairs (empty)\n", grp))
  }

  # Pathway-level table
  lr_pathway <- subsetCommunication(cc, slot.name = "netP")
  if (nrow(lr_pathway) > 0) {
    fwrite(lr_pathway, file.path(out$reports, sprintf("pathway_communication_%s.csv", grp)))
    cat(sprintf("  %s: %d pathway-level entries\n", grp, nrow(lr_pathway)))
  }
}

# ============================================================================
# 6. Collect inference summary (in-memory only; export controlled by save_qc_tables)
# ============================================================================
summary_df <- do.call(rbind, inference_summary)

cat("\n[QC] Inference summary:\n")
print(summary_df, row.names = FALSE)


# ============================================================================
# 7. Provenance
# ============================================================================
if (isTRUE(cfg$output$save_config_copy)) {
  file.copy(config_path, file.path(out$qc, "cellchat_params_snapshot_02.yaml"), overwrite = TRUE)
}

if (isTRUE(cfg$output$save_session_info)) {
  writeLines(capture.output(sessionInfo()),
             file.path(out$qc, "02_run_session_info.txt"))
}

# ── save_cellchat_version：与 Step 01 保持一致 ──
if (isTRUE(cfg$output$save_cellchat_version)) {
  writeLines(
    c(sprintf("CellChat: %s", packageVersion("CellChat")),
      sprintf("R:        %s", R.version.string),
      sprintf("Dataset:  %s", ds_id),
      sprintf("Step:     02_run_cellchat"),
      sprintf("Date:     %s", Sys.time())),
    file.path(out$qc, "versions_02.txt")
  )
}

# ── save_qc_tables：控制 inference_summary 导出 ──
if (isTRUE(cfg$output$save_qc_tables)) {
  fwrite(summary_df, file.path(out$qc, "inference_summary.csv"))
  cat(sprintf("[qc] inference_summary.csv written (%d rows)\n", nrow(summary_df)))
} else {
  cat("[qc] save_qc_tables=false → inference_summary.csv skipped\n")
}


# ============================================================================
# 8. Summary
# ============================================================================
cat("\n")
rule <- strrep("=", 64)
cat(sprintf("+%s+\n", rule))
cat("| CellChat Step 02 Complete\n")
cat(sprintf("+%s+\n", rule))
cat(sprintf("| Dataset:    %s\n", ds_id))
cat(sprintf("| Inference:  type=%s, nboot=%d, thresh=%.3f, pop.size=%s\n",
            inf$type, inf$nboot, inf$thresh, pop_size))
cat(sprintf("+%s+\n", rule))

for (grp in analysis_groups) {
  s <- inference_summary[[grp]]
  cat(sprintf("| %-8s  %5d cells  %2d celltypes  %5d L-R  %3d pathways  %.1f min\n",
              s$group, s$n_cells, s$n_celltypes, s$n_LR_pairs, s$n_pathways, s$elapsed_min))
}

total_min <- sum(summary_df$elapsed_min)
cat(sprintf("+%s+\n", rule))
cat(sprintf("| Total time: %.1f min\n", total_min))
cat(sprintf("| Output:     %s\n", out$base))

# Dynamic next-step guidance
next_script <- file.path(module_root, "scripts", "03_compare_groups.R")
if (file.exists(next_script)) {
  cat("| Next:       Rscript scripts/03_compare_groups.R\n")
} else {
  cat("| Next:       03_compare_groups.R  (not yet created)\n")
}
cat(sprintf("+%s+\n", rule))
