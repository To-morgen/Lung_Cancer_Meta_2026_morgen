#!/usr/bin/env Rscript
# ============================================================================
# 03_recover_scevan_artifacts.R — Priority-based recovery from SCEVAN .RData
#
# Strategy:
#   1. Load LLC_tumor_CNAmtxSubclones.RData
#   2. Extract results.com → tumor_assigned_subclone (~17827 cells)
#   3. Remaining CNA_mtx_relat cells minus ref → tumor_unassigned_subclone (~228)
#   4. All other Seurat cells → non_tumor
#   5. Plausibility guard before writing any output
#
# NEVER uses CNA_mtx_relat alone as tumor source (that's ALL ~60k cells)
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(data.table)
  library(yaml)
  library(ggplot2)
  library(patchwork)
})

module_root <- here::here()
proj_root   <- Sys.getenv("LUNGMETA_ROOT")
if (proj_root == "") proj_root <- normalizePath(file.path(module_root, "..", ".."))

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

path_from_root <- function(path) {
  if (grepl("^/", path)) path else file.path(proj_root, path)
}

cat("\n")
cat("╔══════════════════════════════════════════════════╗\n")
cat("║   03: SCEVAN Artifact Recovery (Priority-based)  ║\n")
cat(sprintf("║   Time: %s            ║\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("╚══════════════════════════════════════════════════╝\n\n")

# ============================================================================
# 1. Config
# ============================================================================

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

seurat_path <- path_from_root(cfg$input$seurat_object)
cluster_col <- cfg$input$cluster_column %||% "seurat_clusters"
ref_clusters <- as.character(cfg$scevan$reference_clusters)

scevan_out <- if (!is.null(cfg$output$scevan)) {
  path_from_root(cfg$output$scevan)
} else {
  file.path(path_from_root(cfg$output$main_results), "scevan")
}
report_out <- if (!is.null(cfg$output$reports)) {
  path_from_root(cfg$output$reports)
} else {
  file.path(path_from_root(cfg$output$main_results), "reports")
}
plot_out <- if (!is.null(cfg$output$plots)) {
  path_from_root(cfg$output$plots)
} else {
  file.path(path_from_root(cfg$output$main_results), "plots")
}
local_out <- Sys.getenv("CNV_LOCAL_SCEVAN_DIR", unset = "")
if (local_out == "") {
  local_out <- if (!is.null(cfg$output$module_scevan)) {
    path_from_root(cfg$output$module_scevan)
  } else {
    file.path(module_root, cfg$output$module_results %||% "results", "scevan")
  }
}
for (d in c(scevan_out, report_out, plot_out, local_out)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat(sprintf("Seurat:     %s\n", basename(seurat_path)))
cat(sprintf("Cluster:    %s\n", cluster_col))
cat(sprintf("Ref clusters: %s\n", paste(ref_clusters, collapse = ", ")))

# ============================================================================
# 2. Load Seurat
# ============================================================================

cat("\nLoading Seurat...\n")
sobj <- readRDS(seurat_path)
cell_barcodes <- colnames(sobj)
n_total <- length(cell_barcodes)
cat(sprintf("  %d cells, %d genes\n", ncol(sobj), nrow(sobj)))

# ============================================================================
# 3. Load .RData artifacts
# ============================================================================

# Priority: Subclones file (has results.com)
scevan_sample <- cfg$scevan$sample_name %||% cfg$dataset_id %||% "LLC_tumor"
subclone_path <- file.path(local_out, sprintf("%s_CNAmtxSubclones.RData", scevan_sample))
base_path     <- file.path(local_out, sprintf("%s_CNAmtx.RData", scevan_sample))
if (!file.exists(subclone_path)) {
  candidates <- list.files(local_out, pattern = "CNAmtxSubclones\\.RData$", full.names = TRUE, recursive = TRUE)
  if (length(candidates) > 0) subclone_path <- candidates[1]
}
if (!file.exists(base_path)) {
  candidates <- list.files(local_out, pattern = "CNAmtx\\.RData$", full.names = TRUE, recursive = TRUE)
  if (length(candidates) > 0) base_path <- candidates[1]
}

if (!file.exists(subclone_path)) {
  stop("CRITICAL: LLC_tumor_CNAmtxSubclones.RData not found at: ", subclone_path)
}

cat(sprintf("\nLoading: %s\n", basename(subclone_path)))
env_sub <- new.env()
load(subclone_path, envir = env_sub)
cat(sprintf("  Objects: %s\n", paste(ls(env_sub), collapse = ", ")))

# Also load base .RData for CNA_mtx_relat (if needed for unassigned)
env_base <- new.env()
if (file.exists(base_path)) {
  load(base_path, envir = env_base)
  cat(sprintf("  Base objects: %s\n", paste(ls(env_base), collapse = ", ")))
}

# ============================================================================
# 4. Extract results.com (PRIMARY tumor source)
# ============================================================================

if (!"results.com" %in% ls(env_sub)) {
  stop("CRITICAL: results.com not found in Subclones .RData!")
}

results_com <- env_sub$results.com

cat(sprintf("\nresults.com: %d rows\n", nrow(results_com)))
cat(sprintf("  Columns: %s\n", paste(colnames(results_com), collapse = ", ")))

# Get barcodes from results.com
rc_barcodes <- colnames(results_com)
# (fallback removed — colnames is the correct source)

# Match to Seurat
rc_in_seurat <- rc_barcodes[rc_barcodes %in% cell_barcodes]
cat(sprintf("  Matched to Seurat: %d / %d\n", length(rc_in_seurat), length(rc_barcodes)))

# These are TUMOR cells with SUBCLONE assignment
tumor_assigned <- rc_in_seurat
n_assigned <- length(tumor_assigned)
cat(sprintf("\n  → tumor_assigned_subclone: %d cells\n", n_assigned))

# ============================================================================
# 5. Identify tumor_unassigned (in CNA_mtx_relat but NOT in results.com, NOT ref)
# ============================================================================

# Get CNA_mtx_relat barcodes (all cells that entered SCEVAN's CNV analysis)
cna_barcodes <- character(0)
for (env in list(env_sub, env_base)) {
  if ("CNA_mtx_relat" %in% ls(env)) {
    cna_barcodes <- colnames(env$CNA_mtx_relat)
    cat(sprintf("\nCNA_mtx_relat: %d columns\n", length(cna_barcodes)))
    break
  }
}

# Reference barcodes (from Seurat cluster metadata)
ref_bc <- cell_barcodes[as.character(sobj@meta.data[[cluster_col]]) %in% ref_clusters]
cat(sprintf("Reference cells (clusters %s): %d\n",
            paste(ref_clusters, collapse = ","), length(ref_bc)))

# SCEVAN classified only results.com colnames as tumor (17,827 of ~18,055).
# The ~228 difference were filtered during subclone refinement → treat as non-tumor.
# All other cells in CNA_mtx_relat were SCEVAN-classified NORMAL — NOT tumor.
tumor_unassigned <- character(0)
n_unassigned <- 0
cat("  → tumor_unassigned: 0 (only results.com colnames are confirmed tumor)\n")

if (length(cna_barcodes) > 0) {
  cna_in_seurat <- cna_barcodes[cna_barcodes %in% cell_barcodes]
  scevan_normal_nonref <- setdiff(cna_in_seurat, c(tumor_assigned, ref_bc))
  cat(sprintf("  → SCEVAN-classified normal (non-ref clusters): %d cells\n",
              length(scevan_normal_nonref)))
}

# Non-tumor = everything else
all_tumor <- union(tumor_assigned, tumor_unassigned)
non_tumor <- setdiff(cell_barcodes, all_tumor)
n_non_tumor <- length(non_tumor)

cat(sprintf("\n=== Label Summary ===\n"))
cat(sprintf("  tumor_assigned_subclone:   %6d  (%.1f%%)\n",
            n_assigned, n_assigned / n_total * 100))
cat(sprintf("  tumor_unassigned_subclone: %6d  (%.1f%%)\n",
            n_unassigned, n_unassigned / n_total * 100))
cat(sprintf("  non_tumor:                 %6d  (%.1f%%)\n",
            n_non_tumor, n_non_tumor / n_total * 100))
cat(sprintf("  TOTAL:                     %6d\n", n_total))

# ============================================================================
# 6. Plausibility Guard (BEFORE writing anything)
# ============================================================================

tumor_frac <- (n_assigned + n_unassigned) / n_total * 100

# Guard 1: tumor fraction > 90% → ABORT (CNA_mtx_relat contamination signature)
if (tumor_frac > 90) {
  stop(sprintf(
    "PLAUSIBILITY FAIL: tumor fraction = %.1f%% (> 90%%). " %+%
    "This looks like CNA_mtx_relat was used instead of results.com. ABORTING.",
    tumor_frac))
}

# Guard 2: ref contamination > 5% → WARNING
ref_in_tumor <- sum(all_tumor %in% ref_bc)
ref_contam <- ref_in_tumor / length(ref_bc) * 100
if (ref_contam > 5) {
  cat(sprintf("⚠️  WARNING: %.1f%% of reference cells labeled as tumor!\n", ref_contam))
} else {
  cat(sprintf("✅ Reference contamination: %.2f%% (acceptable)\n", ref_contam))
}

# Guard 3: tumor fraction < 5% → WARNING
if (tumor_frac < 5) {
  cat(sprintf("⚠️  WARNING: tumor fraction = %.1f%% — unusually low\n", tumor_frac))
}

cat(sprintf("✅ Plausibility check passed: tumor = %.1f%%\n", tumor_frac))

# ============================================================================
# 7. Build label table
# ============================================================================

label_df <- data.table(
  barcode    = cell_barcodes,
  scevan_call = "non_tumor"
)
label_df[barcode %in% tumor_assigned, scevan_call := "tumor_assigned_subclone"]
label_df[barcode %in% tumor_unassigned, scevan_call := "tumor_unassigned_subclone"]

# Add subclone info from results.com (if available)
subcl_cols <- grep("subclone|class", colnames(results_com), ignore.case = TRUE, value = TRUE)
if (length(subcl_cols) > 0) {
  # Build subclone lookup
  subcl_lookup <- data.table(
    barcode = rc_barcodes[rc_barcodes %in% cell_barcodes]
  )
  for (sc in subcl_cols) {
    vals <- results_com[[sc]]
    if (length(vals) == length(rc_barcodes)) {
      subcl_lookup[[sc]] <- vals[rc_barcodes %in% cell_barcodes]
    }
  }
  label_df <- merge(label_df, subcl_lookup, by = "barcode", all.x = TRUE)
}

# Add simplified binary call
label_df[, is_tumor := grepl("^tumor", scevan_call)]

cat(sprintf("\n=== Final Label Distribution ===\n"))
print(label_df[, .N, by = scevan_call])

# ============================================================================
# 8. Add to Seurat + per-cluster/group summaries
# ============================================================================

# Raw labels CSV
fwrite(label_df, file.path(scevan_out, "scevan_recovered_labels_raw.csv"))
fwrite(label_df, file.path(report_out, "scevan_recovered_labels_full.csv"))

# Add cluster + group info
md <- sobj@meta.data
label_df[, cluster := as.character(md[barcode, cluster_col])]
label_df[, group := as.character(md[barcode, "group"])]
label_df[, sample_id := as.character(md[barcode, "sample_id"])]

# Per-cluster
cluster_summary <- label_df[, .(
  n_total  = .N,
  n_tumor  = sum(is_tumor),
  n_assigned = sum(scevan_call == "tumor_assigned_subclone"),
  n_unassigned = sum(scevan_call == "tumor_unassigned_subclone"),
  pct_tumor = round(sum(is_tumor) / .N * 100, 2),
  is_reference = unique(cluster) %in% ref_clusters
), by = cluster][order(as.integer(cluster))]

fwrite(cluster_summary, file.path(report_out, "scevan_recovered_tumor_by_cluster.csv"))
cat("\n=== Tumor by Cluster ===\n")
print(cluster_summary)

# Per-group
group_summary <- label_df[, .(
  n_total  = .N,
  n_tumor  = sum(is_tumor),
  pct_tumor = round(sum(is_tumor) / .N * 100, 2)
), by = group]

fwrite(group_summary, file.path(report_out, "scevan_recovered_tumor_by_group.csv"))
cat("\n=== Tumor by Group ===\n")
print(group_summary)

# Per-sample
sample_summary <- label_df[, .(
  n_total  = .N,
  n_tumor  = sum(is_tumor),
  pct_tumor = round(sum(is_tumor) / .N * 100, 2)
), by = sample_id]

fwrite(sample_summary, file.path(report_out, "scevan_recovered_tumor_by_sample.csv"))
cat("\n=== Tumor by Sample ===\n")
print(sample_summary)

# Reference contamination check
ref_check <- label_df[cluster %in% ref_clusters, .(
  n_total = .N,
  n_tumor = sum(is_tumor),
  pct_tumor = round(sum(is_tumor) / .N * 100, 4)
), by = cluster][order(as.integer(cluster))]

fwrite(ref_check, file.path(report_out, "scevan_recovered_reference_check.csv"))
cat("\n=== Reference Cluster Contamination ===\n")
print(ref_check)

# ============================================================================
# 9. Save labeled Seurat
# ============================================================================

cat("\nAdding labels to Seurat...\n")

# Create metadata columns
meta_add <- data.frame(row.names = label_df$barcode)
meta_add$scevan_call <- label_df$scevan_call
meta_add$scevan_is_tumor <- label_df$is_tumor
meta_add$scevan_label <- ifelse(label_df$is_tumor, "tumor", "non_tumor")
meta_add$scevan_source <- "recovery"
meta_add$scevan_subclone <- NA_character_

# Add subclone columns if present
for (sc in subcl_cols) {
  if (sc %in% colnames(label_df)) {
    meta_add[[paste0("scevan_", sc)]] <- label_df[[sc]]
    if (grepl("subclone", sc, ignore.case = TRUE)) {
      meta_add$scevan_subclone <- as.character(label_df[[sc]])
    }
  }
}

sobj <- AddMetaData(sobj, metadata = meta_add)

out_rds <- file.path(scevan_out, "seurat_with_scevan_recovered.rds")
saveRDS(sobj, out_rds)
cat(sprintf("  Seurat → %s\n", out_rds))

# ============================================================================
# 10. Plots
# ============================================================================

cat("\nGenerating plots...\n")

# Plot 1: UMAP tumor/normal
tryCatch({
  tumor_cols <- c("tumor_assigned_subclone" = "#E41A1C",
                  "tumor_unassigned_subclone" = "#FF7F00",
                  "non_tumor" = "#4DAF4A")

  p1 <- DimPlot(sobj, group.by = "scevan_call", reduction = "umap",
                cols = tumor_cols, pt.size = 0.2) +
    labs(title = sprintf("SCEVAN Recovery: %d tumor / %d total",
                         n_assigned + n_unassigned, n_total))

  p2 <- DimPlot(sobj, group.by = "scevan_call", reduction = "umap",
                split.by = "group", cols = tumor_cols, pt.size = 0.2) +
    labs(title = "By Group")

  pdf(file.path(plot_out, "01_scevan_recovered_umap.pdf"), width = 18, height = 6)
  print(p1 | p2)
  dev.off()

  png(file.path(plot_out, "01_scevan_recovered_umap.png"), width = 1800, height = 600, res = 150)
  print(p1 | p2)
  dev.off()

  cat("  ✅ 01_scevan_recovered_umap\n")
}, error = function(e) cat(sprintf("  ⚠️ UMAP failed: %s\n", e$message)))

# Plot 2: Tumor fraction bar chart
tryCatch({
  p3 <- ggplot(cluster_summary, aes(x = reorder(cluster, -pct_tumor), y = pct_tumor)) +
    geom_col(aes(fill = ifelse(is_reference, "reference", ifelse(pct_tumor > 50, "tumor_dominant", "mixed"))),
             show.legend = TRUE) +
    scale_fill_manual(values = c("reference" = "#377EB8", "tumor_dominant" = "#E41A1C", "mixed" = "#FF7F00"),
                      name = "Type") +
    geom_hline(yintercept = 50, linetype = "dashed", color = "grey40") +
    geom_text(aes(label = sprintf("%.0f%%", pct_tumor)), vjust = -0.3, size = 2.5) +
    labs(title = "Tumor Fraction by Cluster (SCEVAN recovered)",
         x = "Cluster", y = "% Tumor") +
    theme_minimal()

  p4 <- ggplot(group_summary, aes(x = group, y = pct_tumor, fill = group)) +
    geom_col(show.legend = FALSE, alpha = 0.8) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct_tumor, n_tumor)),
              vjust = -0.1, size = 3.5) +
    labs(title = "Tumor Fraction by Group", x = "", y = "% Tumor") +
    theme_minimal()

  pdf(file.path(plot_out, "02_scevan_recovered_tumor_fraction.pdf"), width = 16, height = 6)
  print(p3 | p4)
  dev.off()

  png(file.path(plot_out, "02_scevan_recovered_tumor_fraction.png"), width = 1600, height = 600, res = 150)
  print(p3 | p4)
  dev.off()

  cat("  ✅ 02_scevan_recovered_tumor_fraction\n")
}, error = function(e) cat(sprintf("  ⚠️ Bar chart failed: %s\n", e$message)))

# ============================================================================
# 11. Run Note
# ============================================================================

run_note <- data.table(
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  source_file = basename(subclone_path),
  primary_source = "results.com",
  n_tumor_assigned = n_assigned,
  n_tumor_unassigned = n_unassigned,
  n_non_tumor = n_non_tumor,
  n_total = n_total,
  tumor_pct = round(tumor_frac, 2),
  ref_contamination_pct = round(ref_contam, 4),
  plausibility = "PASS"
)
fwrite(run_note, file.path(report_out, "scevan_recovery_run_note.csv"))

gc()

cat("\n")
cat("╔══════════════════════════════════════════════════════╗\n")
cat("║            Recovery Complete                         ║\n")
cat(sprintf("║  Tumor assigned:   %6d  (%.1f%%)                   ║\n",
            n_assigned, n_assigned / n_total * 100))
cat(sprintf("║  Tumor unassigned: %6d  (%.1f%%)                    ║\n",
            n_unassigned, n_unassigned / n_total * 100))
cat(sprintf("║  Non-tumor:        %6d  (%.1f%%)                   ║\n",
            n_non_tumor, n_non_tumor / n_total * 100))
cat(sprintf("║  Ref contamination: %.2f%%                          ║\n", ref_contam))
cat("╠══════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Seurat:  %s\n", out_rds))
cat(sprintf("║  Reports: %s\n", report_out))
cat(sprintf("║  Plots:   %s\n", plot_out))
cat("╚══════════════════════════════════════════════════════╝\n")
