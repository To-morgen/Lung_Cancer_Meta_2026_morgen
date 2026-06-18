# ============================================================================
# 01_prepare_cellchat.R — Prepare per-group CellChat objects
#
# Module:  modules/cellchat/  (stateless — results go to project-level)
# Input:   seurat_annotated_final.rds (Source of Truth)
# Output:  results/scrna/09_cellchat/{dataset_id}/objects/  +  qc/
#
# Usage:   cd modules/cellchat && Rscript scripts/01_prepare_cellchat.R
# Override dataset: DATASET_ID=gse253718_luad Rscript scripts/01_prepare_cellchat.R
# ============================================================================

cat("
╔══════════════════════════════════════════════════════════════╗
║        CellChat Step 01: Prepare Per-Group Objects          ║
╚══════════════════════════════════════════════════════════════╝
\n")

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(Seurat)
  library(CellChat)
  library(data.table)
})

# ============================================================================
# 0. Path resolution
# ============================================================================
module_root  <- here::here()
project_root <- normalizePath(file.path(module_root, "..", ".."))

config_path <- file.path(module_root, "configs", "cellchat_params.yaml")
if (!file.exists(config_path)) stop("Config not found: ", config_path)
cfg <- yaml::read_yaml(config_path)
cat(sprintf("[config] Loaded: %s\n", config_path))
cat(sprintf("[paths]  module_root:  %s\n", module_root))
cat(sprintf("[paths]  project_root: %s\n", project_root))

# ============================================================================
# Helper: resolve_dataset_id
# ============================================================================
resolve_dataset_id <- function(cfg) {
  ds_id <- Sys.getenv("DATASET_ID", "")
  if (ds_id == "" && !is.null(cfg$dataset$id)) ds_id <- cfg$dataset$id
  if (ds_id == "") stop("[CONTRACT] No dataset_id: set DATASET_ID env var or config dataset.id")
  cat(sprintf("[dataset] id: %s\n", ds_id))
  ds_id
}

# ============================================================================
# Helper: resolve_output_dirs  (dataset-first, {dataset_id} template)
# ============================================================================
resolve_output_dirs <- function(cfg, project_root, dataset_id) {
  base_template <- cfg$output$base_dir
  base_resolved <- gsub("\\{dataset_id\\}", dataset_id, base_template)
  base <- file.path(project_root, base_resolved)

  dirs <- list(
    base    = base,
    objects = file.path(base, "objects"),
    figures = file.path(base, "figures"),
    qc      = file.path(base, "qc"),
    reports = file.path(base, "reports")
  )
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("[dirs] canonical output: %s\n", base))
  dirs
}

# ============================================================================
# Helper: resolve_design
# ============================================================================
resolve_design <- function(cfg, project_root, sobj = NULL) {
  src <- cfg$design$source
  cat(sprintf("[design] source: %s\n", src))

  if (src == "de_contrasts") {
    de_path_raw <- Sys.getenv("DE_CONTRASTS_CONFIG", "")
    if (de_path_raw == "") de_path_raw <- cfg$design$contrasts_config
    if (is.null(de_path_raw) || de_path_raw == "") {
      stop("[CONTRACT] No contrasts config: set DE_CONTRASTS_CONFIG or design.contrasts_config")
    }
    de_path <- if (grepl("^/", de_path_raw)) de_path_raw else file.path(project_root, de_path_raw)
    if (!file.exists(de_path)) stop("[CONTRACT] Contrasts config not found: ", de_path)
    de_cfg <- yaml::read_yaml(de_path)

    all_groups <- unique(unlist(lapply(de_cfg$contrasts, function(c) {
      c(c$numerator, c$denominator)
    })))
    reference  <- de_cfg$reference_group
    comparisons <- lapply(names(de_cfg$contrasts), function(nm) {
      list(
        name   = nm,
        groups = c(de_cfg$contrasts[[nm]]$numerator,
                   de_cfg$contrasts[[nm]]$denominator)
      )
    })
    cat(sprintf("[design] Loaded: %s\n", de_path))

  } else if (src == "auto_pairwise") {
    if (is.null(sobj)) stop("[CONTRACT] auto_pairwise requires sobj argument")
    grps <- sort(unique(sobj@meta.data[[cfg$input$group_col]]))
    reference <- grps[1]
    pairs <- combn(grps, 2, simplify = FALSE)
    comparisons <- lapply(pairs, function(p) {
      if (p[2] == reference) p <- rev(p)
      list(name = paste(p[1], "vs", p[2], sep = "_"), groups = p)
    })
    all_groups <- grps
    cat("[design] Auto-generated from Seurat metadata\n")

  } else if (src == "inline") {
    all_groups  <- cfg$design$inline_groups
    reference   <- cfg$design$inline_reference
    comparisons <- cfg$design$inline_comparisons
    if (length(all_groups) == 0) stop("[CONTRACT] inline design but inline_groups is empty")
    cat("[design] Using inline definitions\n")

  } else {
    stop("[CONTRACT] Unknown design source: ", src)
  }

  cat(sprintf("[design] Groups: %s\n", paste(all_groups, collapse = ", ")))
  cat(sprintf("[design] Reference: %s\n", reference))
  cat(sprintf("[design] Comparisons: %s\n",
              paste(sapply(comparisons, `[[`, "name"), collapse = ", ")))

  list(
    all_groups      = all_groups,
    reference_group = reference,
    comparisons     = comparisons
  )
}

# ============================================================================
# Helper: validate_metadata_cols
# ============================================================================
validate_metadata_cols <- function(sobj, cols, context = "Seurat") {
  available <- colnames(sobj@meta.data)
  for (col in cols) {
    if (!col %in% available) {
      stop(sprintf(
        "[CONTRACT] Column '%s' not found in %s metadata.\n  Available: %s",
        col, context,
        paste(head(sort(available), 30), collapse = ", ")
      ))
    }
  }
  cat(sprintf("[validate] Required columns confirmed: %s\n", paste(cols, collapse = ", ")))
}

# ============================================================================
# 1. Resolve dataset & output paths
# ============================================================================
dataset_id <- resolve_dataset_id(cfg)
out        <- resolve_output_dirs(cfg, project_root, dataset_id)

# ============================================================================
# 2. Load Seurat object
# ============================================================================
seurat_path_template <- cfg$input$seurat_path
seurat_path_resolved <- gsub("\\{dataset_id\\}", dataset_id, seurat_path_template)
seurat_path <- file.path(project_root, seurat_path_resolved)
cat(sprintf("[load] Resolved path: %s\n", seurat_path))

if (!file.exists(seurat_path)) {
  stop(sprintf("[CONTRACT] Seurat object not found: %s\n  Check input.seurat_path in config and dataset_id '%s'",
               seurat_path, dataset_id))
}

sobj <- readRDS(seurat_path)
cat(sprintf("[load] %d cells × %d genes\n", ncol(sobj), nrow(sobj)))

celltype_col <- cfg$input$celltype_col
group_col    <- cfg$input$group_col

# ── Input contract validation ──
validate_metadata_cols(sobj, c(celltype_col, group_col))
if (!is.null(cfg$filtering$exclude_L1)) {
  validate_metadata_cols(sobj, "celltype_L1")
}

cat(sprintf("[load] Groups: %s\n",
            paste(sort(unique(sobj@meta.data[[group_col]])), collapse = ", ")))
cat(sprintf("[load] Celltypes (%s): %d unique\n",
            celltype_col, length(unique(sobj@meta.data[[celltype_col]]))))

# ============================================================================
# 3. Resolve experimental design
# ============================================================================
design <- resolve_design(cfg, project_root, sobj)

# Validate groups exist in Seurat
groups_in_data <- unique(sobj@meta.data[[group_col]])
missing_groups <- setdiff(design$all_groups, groups_in_data)
if (length(missing_groups) > 0) {
  stop(sprintf("[CONTRACT] Groups in design but not in Seurat: %s\n  Available: %s",
               paste(missing_groups, collapse = ", "),
               paste(groups_in_data, collapse = ", ")))
}

# ============================================================================
# 4. Filter cells
# ============================================================================
n_before <- ncol(sobj)

# 4a. Exclude by L1
if (!is.null(cfg$filtering$exclude_L1) && "celltype_L1" %in% colnames(sobj@meta.data)) {
  exclude_L1 <- cfg$filtering$exclude_L1
  mask_L1 <- sobj@meta.data$celltype_L1 %in% exclude_L1
  n_exc <- sum(mask_L1)
  if (n_exc > 0) {
    cat(sprintf("[filter] Excluding L1 = {%s}: %d cells\n",
                paste(exclude_L1, collapse = ", "), n_exc))
    sobj <- sobj[, !mask_L1]
  }
}

# 4b. Exclude by L2
if (!is.null(cfg$filtering$exclude_L2)) {
  exclude_L2 <- cfg$filtering$exclude_L2
  mask_L2 <- sobj@meta.data[[celltype_col]] %in% exclude_L2
  n_exc <- sum(mask_L2)
  if (n_exc > 0) {
    cat(sprintf("[filter] Excluding L2 = {%s}: %d cells\n",
                paste(exclude_L2, collapse = ", "), n_exc))
    sobj <- sobj[, !mask_L2]
  }
}

cat(sprintf("[filter] %d → %d cells after exclusion\n", n_before, ncol(sobj)))

# ── Guard: post-exclusion check ──
if (ncol(sobj) == 0) {
  stop("[FATAL] 0 cells remain after L1/L2 exclusion. Check filtering.exclude_L1 and filtering.exclude_L2.")
}

# Clean factor levels
sobj@meta.data[[celltype_col]] <- droplevels(factor(sobj@meta.data[[celltype_col]]))
sobj@meta.data[[group_col]]    <- droplevels(factor(sobj@meta.data[[group_col]]))

# ============================================================================
# 5. Per-group cell count QC table
# ============================================================================
meta <- sobj@meta.data
count_table <- as.data.frame.matrix(table(meta[[celltype_col]], meta[[group_col]]))
count_table <- cbind(celltype = rownames(count_table), count_table)
rownames(count_table) <- NULL

col_order <- c("celltype", intersect(design$all_groups, colnames(count_table)))
extra_cols <- setdiff(colnames(count_table), col_order)
if (length(extra_cols) > 0) col_order <- c(col_order, extra_cols)
count_table <- count_table[, col_order, drop = FALSE]

cat("\n[QC] Per-group cell counts:\n")
print(count_table, row.names = FALSE)

fwrite(count_table, file.path(out$qc, "celltype_group_counts.csv"))

# ============================================================================
# 6. Load CellChat database
# ============================================================================
db_species <- cfg$database$species

CellChatDB <- switch(db_species,
  "mouse" = CellChatDB.mouse,
  "human" = CellChatDB.human,
  stop(sprintf("[CONTRACT] Unsupported species '%s'. Must be 'mouse' or 'human'.", db_species))
)
cat(sprintf("[db] Species: %s\n", db_species))

if (!isTRUE(cfg$database$run_full_database) && !is.null(cfg$database$categories)) {
  cat(sprintf("[db] Subsetting CellChatDB to: %s\n",
              paste(cfg$database$categories, collapse = ", ")))
  CellChatDB_use <- subsetDB(CellChatDB, search = cfg$database$categories)
} else {
  cat("[db] Using full CellChatDB\n")
  CellChatDB_use <- CellChatDB
}
cat(sprintf("[db] %d interactions in database subset\n", nrow(CellChatDB_use$interaction)))

if (isTRUE(cfg$output$save_db_subset)) {
  fwrite(CellChatDB_use$interaction, file.path(out$qc, "cellchatdb_interactions_used.csv"))
}

# ============================================================================
# 7. Determine per-group celltype availability (pre-scan)
# ============================================================================
min_cells  <- cfg$filtering$min_cells_per_celltype
all_groups <- design$all_groups
drop_log   <- list()

# 7a. First pass: which celltypes pass min_cells in each group?
ct_pass_per_group <- list()

for (grp in all_groups) {
  cells_grp <- rownames(meta[meta[[group_col]] == grp, ])

  if (length(cells_grp) == 0) {
    stop(sprintf("[FATAL] Group '%s' has 0 cells after exclusion filtering.", grp))
  }

  ct_counts <- table(meta[cells_grp, celltype_col])
  ct_pass   <- names(ct_counts[ct_counts >= min_cells])
  ct_drop   <- names(ct_counts[ct_counts < min_cells])

  cat(sprintf("\n[pre-scan] %s: %d cells, %d/%d celltypes pass min_cells=%d\n",
              grp, length(cells_grp), length(ct_pass), length(ct_counts), min_cells))

  if (length(ct_drop) > 0) {
    for (ct in ct_drop) {
      cat(sprintf("  drop: %s (%d cells)\n", ct, as.integer(ct_counts[ct])))
      drop_log[[length(drop_log) + 1]] <- data.frame(
        group    = grp,
        celltype = ct,
        n_cells  = as.integer(ct_counts[ct]),
        reason   = sprintf("below min_cells threshold (%d)", min_cells),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(ct_pass) == 0) {
    stop(sprintf(
      "[FATAL] Group '%s' has 0 celltypes with >= %d cells.\n  Counts: %s",
      grp, min_cells,
      paste(sprintf("%s=%d", names(ct_counts), ct_counts), collapse = ", ")
    ))
  }

  ct_pass_per_group[[grp]] <- ct_pass
}

# 7b. Compute common celltypes
common_all <- Reduce(intersect, ct_pass_per_group)
cat(sprintf("\n[consistency] Celltypes common to ALL groups (%d): %s\n",
            length(common_all), paste(common_all, collapse = ", ")))

if (length(common_all) == 0) {
  stop("[FATAL] 0 celltypes are common across all groups. Cannot proceed.")
}

# 7c. Per-comparison consistency log
cat("\n[consistency] Per-comparison breakdown:\n")
consistency_log <- list()

for (comp in design$comparisons) {
  g1 <- comp$groups[1]
  g2 <- comp$groups[2]
  cts_g1  <- ct_pass_per_group[[g1]]
  cts_g2  <- ct_pass_per_group[[g2]]
  common  <- intersect(cts_g1, cts_g2)
  only_g1 <- setdiff(cts_g1, cts_g2)
  only_g2 <- setdiff(cts_g2, cts_g1)

  cat(sprintf("  %s: %d common | %d only-%s | %d only-%s\n",
              comp$name, length(common), length(only_g1), g1, length(only_g2), g2))
  if (length(only_g1) > 0) cat(sprintf("    only-%s: %s\n", g1, paste(only_g1, collapse = ", ")))
  if (length(only_g2) > 0) cat(sprintf("    only-%s: %s\n", g2, paste(only_g2, collapse = ", ")))

  consistency_log[[comp$name]] <- data.frame(
    comparison = comp$name,
    celltype   = c(common, only_g1, only_g2),
    status     = c(rep("common", length(common)),
                   rep(paste0("only_", g1), length(only_g1)),
                   rep(paste0("only_", g2), length(only_g2))),
    stringsAsFactors = FALSE
  )
}

consistency_df <- do.call(rbind, consistency_log)
fwrite(consistency_df, file.path(out$qc, "celltype_consistency_per_comparison.csv"))

# 7d. Determine final celltype set
if (isTRUE(cfg$filtering$drop_if_absent_in_any_group)) {
  final_celltypes <- common_all
  cat(sprintf("\n[filter] drop_if_absent_in_any_group=TRUE → using %d common celltypes\n",
              length(final_celltypes)))

  # Log dropped celltypes per group
  for (grp in all_groups) {
    to_drop <- setdiff(ct_pass_per_group[[grp]], final_celltypes)
    if (length(to_drop) > 0) {
      cat(sprintf("  %s: dropping %d group-specific celltypes: %s\n",
                  grp, length(to_drop), paste(to_drop, collapse = ", ")))
      for (ct in to_drop) {
        cells_grp <- rownames(meta[meta[[group_col]] == grp, ])
        n <- sum(meta[cells_grp, celltype_col] == ct)
        drop_log[[length(drop_log) + 1]] <- data.frame(
          group    = grp,
          celltype = ct,
          n_cells  = as.integer(n),
          reason   = "absent in at least one other group (drop_if_absent_in_any_group)",
          stringsAsFactors = FALSE
        )
      }
    }
  }
} else {
  # Each group keeps its own passing celltypes
  final_celltypes <- NULL  # NULL = use per-group ct_pass
  cat("\n[filter] drop_if_absent_in_any_group=FALSE → each group keeps its own celltypes\n")
}

# Save per-comparison common celltypes for 03_compare_groups.R
comparison_common <- list()
for (comp in design$comparisons) {
  g1_cts <- if (!is.null(final_celltypes)) final_celltypes else ct_pass_per_group[[comp$groups[1]]]
  g2_cts <- if (!is.null(final_celltypes)) final_celltypes else ct_pass_per_group[[comp$groups[2]]]
  comparison_common[[comp$name]] <- intersect(g1_cts, g2_cts)
}
saveRDS(comparison_common, file.path(out$objects, "comparison_common_celltypes.rds"))

# ============================================================================
# 8. Create CellChat objects (Seurat-level filtering → createCellChat)
# ============================================================================
cellchat_list <- list()

for (grp in all_groups) {
  cat(sprintf("\n━━━ Group: %s ━━━\n", grp))

  # Determine which celltypes this group gets
  grp_celltypes <- if (!is.null(final_celltypes)) final_celltypes else ct_pass_per_group[[grp]]

  # Subset Seurat: group + passing celltypes
  cells_grp <- rownames(meta[meta[[group_col]] == grp &
                               meta[[celltype_col]] %in% grp_celltypes, ])

  if (length(cells_grp) == 0) {
    stop(sprintf("[FATAL] Group '%s' has 0 cells after celltype filtering.", grp))
  }

  sobj_grp <- sobj[, cells_grp]
  sobj_grp@meta.data[[celltype_col]] <- droplevels(factor(sobj_grp@meta.data[[celltype_col]]))


      # Map sample_id → samples (CellChat v2 convention)
  if ("sample_id" %in% colnames(sobj_grp@meta.data) &&
      !"samples" %in% colnames(sobj_grp@meta.data)) {
    sobj_grp@meta.data$samples <- sobj_grp@meta.data$sample_id
  }

  cat(sprintf("  Cells: %d | Celltypes: %d\n",
              ncol(sobj_grp), length(levels(sobj_grp@meta.data[[celltype_col]]))))

  # Create CellChat object (v2 API)
  cellchat <- createCellChat(
    object   = sobj_grp,
    group.by = celltype_col,
    assay    = "RNA"
  )
  cellchat@DB <- CellChatDB_use

  cellchat_list[[grp]] <- cellchat
  cat(sprintf("  CellChat object created: %d cells, %d celltypes\n",
              length(cellchat@idents), length(levels(cellchat@idents))))
}


# ============================================================================
# 9. Save objects
# ============================================================================
cat("\n[save] Saving CellChat objects...\n")

for (grp in names(cellchat_list)) {
  fpath <- file.path(out$objects, sprintf("cellchat_%s.rds", grp))
  saveRDS(cellchat_list[[grp]], fpath)
  cat(sprintf("  %s\n", fpath))
}

list_path <- file.path(out$objects, "cellchat_list.rds")
saveRDS(cellchat_list, list_path)
cat(sprintf("  %s\n", list_path))

design_path <- file.path(out$objects, "design.rds")
saveRDS(design, design_path)
cat(sprintf("  %s\n", design_path))

# ============================================================================
# 10. Save drop log
# ============================================================================
if (length(drop_log) > 0) {
  drop_df <- do.call(rbind, drop_log)
  fwrite(drop_df, file.path(out$qc, "dropped_celltypes.csv"))
  cat(sprintf("[qc] %d celltype×group entries dropped → dropped_celltypes.csv\n", nrow(drop_df)))
} else {
  cat("[qc] No celltypes dropped\n")
}

# ============================================================================
# 11. Provenance
# ============================================================================
if (isTRUE(cfg$output$save_config_copy)) {
  file.copy(config_path, file.path(out$qc, "cellchat_params_snapshot.yaml"), overwrite = TRUE)
}

if (isTRUE(cfg$output$save_session_info)) {
  writeLines(capture.output(sessionInfo()),
             file.path(out$qc, "01_prepare_session_info.txt"))
}

if (isTRUE(cfg$output$save_cellchat_version)) {
  writeLines(
    c(sprintf("CellChat: %s", packageVersion("CellChat")),
      sprintf("Seurat:   %s", packageVersion("Seurat")),
      sprintf("R:        %s", R.version.string),
      sprintf("Dataset:  %s", dataset_id),
      sprintf("Date:     %s", Sys.time())),
    file.path(out$qc, "versions.txt")
  )
}

# ============================================================================
# 12. Summary
# ============================================================================
cat("\n")
rule <- strrep("=", 64)
cat(sprintf("+%s+\n", rule))
cat(sprintf("| CellChat Step 01 Complete\n"))
cat(sprintf("+%s+\n", rule))
cat(sprintf("| Dataset:     %s\n", dataset_id))
cat(sprintf("| Design src:  %s\n", cfg$design$source))
cat(sprintf("| Species:     %s\n", db_species))
cat(sprintf("| DB:          %d interactions (%s)\n",
            nrow(CellChatDB_use$interaction),
            if (isTRUE(cfg$database$run_full_database)) "full" else "subset"))
cat(sprintf("+%s+\n", rule))
for (grp in names(cellchat_list)) {
  cc <- cellchat_list[[grp]]
  cat(sprintf("| %-8s  %5d cells  %2d celltypes\n",
              grp, length(cc@idents), length(levels(cc@idents))))
}
cat(sprintf("+%s+\n", rule))
cat(sprintf("| Common celltypes (all groups): %d\n", length(common_all)))
cat(sprintf("| Comparisons: %s\n",
            paste(sapply(design$comparisons, `[[`, "name"), collapse = ", ")))
cat(sprintf("+%s+\n", rule))
cat(sprintf("| Output: %s\n", out$base))
cat(sprintf("| Next:   Rscript scripts/02_run_cellchat.R\n"))
cat(sprintf("+%s+\n", rule))
