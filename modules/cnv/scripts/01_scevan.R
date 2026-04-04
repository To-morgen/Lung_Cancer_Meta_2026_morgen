#!/usr/bin/env Rscript
# ============================================================================
# 01_scevan.R — SCEVAN CNV inference: tumor/normal classification
#
# Environment: modules/cnv/ (独立 renv)
# Input:  主项目 annotated Seurat object (via LUNGMETA_ROOT)
# Output: tumor/normal labels, CNV heatmap, per-cluster/group summary
#
# Usage:
#   cd ~/biohub/projects/Lung_Cancer_Meta_2026_morgen/modules/cnv
#   Rscript scripts/01_scevan.R
#
# Key improvements over old script:
#   - NO setwd()           → all paths are absolute via file.path()
#   - NO as.matrix()       → SCEVAN accepts sparse; saves ~10GB RAM
#   - NO future setup      → SCEVAN uses parallel::mclapply internally
#   - Config-driven        → reference clusters, cores, organism from YAML
#   - Structured output    → semantic file names, CSV summaries
# ============================================================================

suppressPackageStartupMessages({
  library(here)        # anchored to modules/cnv/ (via cnv.Rproj)
  library(Seurat)
  library(SCEVAN)
  library(data.table)
  library(yaml)
  library(ggplot2)
  library(patchwork)
})

# ============================================================================
# 1. Paths & Config
# ============================================================================

proj_root <- Sys.getenv("LUNGMETA_ROOT")
if (proj_root == "") {
  # Fallback: infer from module location
  proj_root <- normalizePath(file.path(here::here(), "..", ".."))
  cat(sprintf("LUNGMETA_ROOT not set, inferred: %s\n", proj_root))
}

module_root <- here::here()  # modules/cnv/

cat("\n")
cat("╔══════════════════════════════════════════════════╗\n")
cat("║   CNV Module: SCEVAN Tumor/Normal Inference      ║\n")
cat(sprintf("║   Project:  %-36s║\n", basename(proj_root)))
cat(sprintf("║   Time:     %s            ║\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("╚══════════════════════════════════════════════════╝\n\n")

# Load config
cfg_path <- file.path(module_root, "configs", "cnv_params.yaml")
if (!file.exists(cfg_path)) stop("Config not found: ", cfg_path)
cfg <- yaml::read_yaml(cfg_path)

cat(sprintf("Species:    %s\n", cfg$species))
cat(sprintf("Organism:   %s\n", cfg$scevan$organism))
cat(sprintf("Cores:      %d\n", cfg$scevan$par_cores))
cat(sprintf("Subclones:  %s\n", ifelse(isTRUE(cfg$scevan$subclones), "yes", "no")))

# Output directories (all absolute — no setwd!)
local_out  <- file.path(module_root, cfg$output$module_results, "scevan")
main_out   <- file.path(proj_root, cfg$output$main_results, "scevan")
plot_out   <- file.path(proj_root, cfg$output$main_results, "plots")
report_out <- file.path(proj_root, cfg$output$main_results, "reports")

for (d in c(local_out, main_out, plot_out, report_out)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================================
# 2. Load Seurat Object
# ============================================================================

input_path <- file.path(proj_root, cfg$input$seurat_object)
if (!file.exists(input_path)) stop("Input not found: ", input_path)

cat(sprintf("\nLoading: %s\n", basename(input_path)))
sobj <- readRDS(input_path)

# Ensure layers are joined (Seurat v5 compatibility)
#sobj <- JoinLayers(sobj)
tryCatch(
  sobj <- JoinLayers(sobj),
  error = function(e) cat(sprintf("  JoinLayers skipped (v4 object): %s\n", e$message))
)


cat(sprintf("  Cells:    %d\n", ncol(sobj)))
cat(sprintf("  Genes:    %d\n", nrow(sobj)))
cat(sprintf("  Clusters: %s\n", paste(sort(unique(sobj$seurat_clusters)), collapse = ", ")))
cat("\n  Cells per sample:\n")
print(table(sobj$sample_id))
cat("\n  Cells per group:\n")
print(table(sobj$group))

# ============================================================================
# 3. Optional Downsample
# ============================================================================

ds_n <- cfg$scevan$downsample
if (!is.null(ds_n) && is.numeric(ds_n) && ds_n < ncol(sobj)) {
  cat(sprintf("\n⚠️  Downsampling: %d → %d cells\n", ncol(sobj), ds_n))
  set.seed(42)
  keep_cells <- sample(colnames(sobj), ds_n)
  sobj <- subset(sobj, cells = keep_cells)
}

# ============================================================================
# 4. Extract Count Matrix (SPARSE — do NOT convert to dense!)
# ============================================================================

cat("\nExtracting count matrix (sparse)...\n")
counts <- GetAssayData(sobj, assay = "RNA", layer = "counts")

# ⚠️ OLD SCRIPT BUG: as.matrix() converts sparse → dense
#    For 50k cells × 30k genes: sparse ~500MB → dense ~12GB
#    SCEVAN handles sparse internally. DO NOT do as.matrix()!

cat(sprintf("  Matrix:   %d genes × %d cells\n", nrow(counts), ncol(counts)))
cat(sprintf("  Class:    %s\n", class(counts)[1]))
cat(sprintf("  Sparse:   %s\n", inherits(counts, "sparseMatrix")))

# Quick sanity check
cat(sprintf("  Total UMI: %.2e\n", sum(counts)))
cat(sprintf("  Median genes/cell: %.0f\n", median(colSums(counts > 0))))

# ============================================================================
# 5. Identify Reference (Normal) Cells
# ============================================================================

ref_config <- cfg$scevan$reference_clusters

if (is.null(ref_config) || (is.character(ref_config) && ref_config == "auto")) {
  cat("\nReference mode: AUTO (SCEVAN will detect normal cells)\n")
  ref_barcodes <- NULL
} else {
  ref_clusters <- as.character(ref_config)
  ref_barcodes <- colnames(sobj)[sobj$seurat_clusters %in% ref_clusters]

  cat(sprintf("\nReference mode: MANUAL\n"))
  cat(sprintf("  Clusters:   %s\n", paste(ref_clusters, collapse = ", ")))
  cat(sprintf("  Ref cells:  %d / %d (%.1f%%)\n",
              length(ref_barcodes), ncol(sobj),
              length(ref_barcodes) / ncol(sobj) * 100))

  if (length(ref_barcodes) < 100) {
    cat("  ⚠️  Warning: < 100 reference cells — consider adding clusters or using 'auto'\n")
  }

  # Log which cell types are in reference
  if ("singler_mouse_cluster" %in% colnames(sobj@meta.data)) {
    ref_md <- sobj@meta.data[sobj$seurat_clusters %in% ref_clusters, ]
    cat("  Reference cell types (SingleR):\n")
    print(table(ref_md$singler_mouse_cluster))
  }
}

# ============================================================================
# 6. Run SCEVAN
# ============================================================================

cat("\n")
cat("══════════════════════════════════════════════════\n")
cat("   Running SCEVAN pipelineCNA\n")
cat(sprintf("   Start: %s\n", format(Sys.time(), "%H:%M:%S")))
cat("══════════════════════════════════════════════════\n\n")

t0 <- Sys.time()

# SCEVAN writes output files to working directory
# Temporarily switch to local output dir, then switch back
old_wd <- getwd()
setwd(local_out)  # ONLY acceptable setwd: for SCEVAN's internal file output

tryCatch({
  if (is.null(ref_barcodes)) {
    # Auto mode
    results <- SCEVAN::pipelineCNA(
      count_mtx  = counts,
      sample     = "LLC_tumor",
      organism   = cfg$scevan$organism,
      par_cores  = cfg$scevan$par_cores,
      SUBCLONES  = isTRUE(cfg$scevan$subclones)
    )
  } else {
    # Manual reference mode
    results <- SCEVAN::pipelineCNA(
      count_mtx  = counts,
      sample     = "LLC_tumor",
      organism   = cfg$scevan$organism,
      par_cores  = cfg$scevan$par_cores,
      norm_cell  = ref_barcodes,
      SUBCLONES  = isTRUE(cfg$scevan$subclones)
    )
  }
}, error = function(e) {
  setwd(old_wd)
  stop("SCEVAN failed: ", e$message)
})

setwd(old_wd)

elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
cat(sprintf("\n✅ SCEVAN completed in %.1f min\n", elapsed))

# ============================================================================
# 7. Parse Results
# ============================================================================

cat("\n========== Parsing SCEVAN output ==========\n")

# SCEVAN returns a data.frame with columns typically including "class"
# But format varies by version — handle robustly

if (is.data.frame(results)) {
  scevan_df <- results
  if (!"barcode" %in% colnames(scevan_df)) {
    scevan_df$barcode <- rownames(scevan_df)
  }
  cat(sprintf("Direct result: %d rows, columns: %s\n",
              nrow(scevan_df), paste(colnames(scevan_df), collapse = ", ")))
} else {
  # Try to find output CSV from SCEVAN's file output
  out_csv <- list.files(local_out, pattern = "LLC_tumor.*\\.csv$",
                        full.names = TRUE, recursive = TRUE)
  if (length(out_csv) > 0) {
    scevan_df <- fread(out_csv[1])
    cat(sprintf("Read from file: %s (%d rows)\n", basename(out_csv[1]), nrow(scevan_df)))
  } else {
    cat("⚠️  Cannot parse SCEVAN output. Saving raw object for manual inspection.\n")
    saveRDS(results, file.path(local_out, "scevan_raw_unparsed.rds"))
    setwd(old_wd)
    stop("SCEVAN output format not recognized. Check scevan_raw_unparsed.rds")
  }
}

# Identify the class column
class_col <- intersect(c("class", "scevan_class", "tumor_normal"), colnames(scevan_df))[1]
if (is.na(class_col)) {
  cat("Columns available: ", paste(colnames(scevan_df), collapse = ", "), "\n")
  cat("⚠️  No recognized class column found\n")
  class_col <- NULL
}

if (!is.null(class_col)) {
  cat(sprintf("\nClass column: '%s'\n", class_col))
  cat("Distribution:\n")
  print(table(scevan_df[[class_col]], useNA = "ifany"))
}

# ============================================================================
# 8. Add Labels to Seurat & Save
# ============================================================================

cat("\n========== Saving outputs ==========\n")

# 8a. Full SCEVAN result table
fwrite(scevan_df, file.path(local_out, "scevan_full_results.csv"))
fwrite(scevan_df, file.path(main_out, "scevan_full_results.csv"))

# 8b. Raw SCEVAN output object (for debugging)
saveRDS(results, file.path(local_out, "scevan_raw_output.rds"))

# 8c. Simplified tumor/normal labels
if (!is.null(class_col)) {
  bc_col <- if ("barcode" %in% colnames(scevan_df)) "barcode" else {
    if (!is.null(rownames(scevan_df)) && !identical(rownames(scevan_df), as.character(1:nrow(scevan_df)))) {
      scevan_df$barcode <- rownames(scevan_df)
      "barcode"
    } else {
      colnames(scevan_df)[1]
    }
  }

  labels <- data.frame(
    barcode     = scevan_df[[bc_col]],
    scevan_call = scevan_df[[class_col]],
    stringsAsFactors = FALSE
  )

  # Add subclone info if exists
  subcl_col <- intersect(c("subclone", "Subclone", "SUBCLONE"), colnames(scevan_df))
  if (length(subcl_col) > 0) {
    labels$subclone <- scevan_df[[subcl_col[1]]]
  }

  fwrite(labels, file.path(main_out, "scevan_tumor_labels.csv"))
  fwrite(labels, file.path(report_out, "scevan_tumor_labels.csv"))

  # 8d. Add to Seurat metadata
  sobj_labeled <- AddMetaData(sobj, metadata = results)
  saveRDS(sobj_labeled, file.path(main_out, "seurat_with_cnv.rds"))
  cat(sprintf("  Seurat + CNV labels → %s\n", file.path(main_out, "seurat_with_cnv.rds")))

  # ============================================================================
  # 9. Summary Tables
  # ============================================================================

  cat("\n========== Summary ==========\n")

  # Match labels to cluster/group metadata
  matched <- merge(
    labels,
    data.frame(
      barcode = colnames(sobj),
      cluster = as.character(sobj$seurat_clusters),
      group   = sobj$group,
      sample  = sobj$sample_id,
      stringsAsFactors = FALSE
    ),
    by = "barcode", all.x = TRUE
  )
  matched_dt <- as.data.table(matched)

  # Per-cluster summary
  cluster_summary <- matched_dt[, .(
    n_total  = .N,
    n_tumor  = sum(scevan_call == "tumor", na.rm = TRUE),
    n_normal = sum(scevan_call == "normal", na.rm = TRUE),
    n_other  = sum(!scevan_call %in% c("tumor", "normal"), na.rm = TRUE),
    pct_tumor = round(sum(scevan_call == "tumor", na.rm = TRUE) / .N * 100, 1)
  ), by = cluster][order(as.integer(cluster))]

  fwrite(cluster_summary, file.path(report_out, "scevan_tumor_by_cluster.csv"))
  cat("\n  Tumor fraction by CLUSTER:\n")
  print(cluster_summary)

  # Per-group summary
  group_summary <- matched_dt[, .(
    n_total  = .N,
    n_tumor  = sum(scevan_call == "tumor", na.rm = TRUE),
    n_normal = sum(scevan_call == "normal", na.rm = TRUE),
    pct_tumor = round(sum(scevan_call == "tumor", na.rm = TRUE) / .N * 100, 1)
  ), by = group]

  fwrite(group_summary, file.path(report_out, "scevan_tumor_by_group.csv"))
  cat("\n  Tumor fraction by GROUP:\n")
  print(group_summary)

  # Per-sample summary
  sample_summary <- matched_dt[, .(
    n_total  = .N,
    n_tumor  = sum(scevan_call == "tumor", na.rm = TRUE),
    pct_tumor = round(sum(scevan_call == "tumor", na.rm = TRUE) / .N * 100, 1)
  ), by = sample]

  fwrite(sample_summary, file.path(report_out, "scevan_tumor_by_sample.csv"))
  cat("\n  Tumor fraction by SAMPLE:\n")
  print(sample_summary)

  # ============================================================================
  # 10. Plots
  # ============================================================================

  cat("\n========== Generating plots ==========\n")

  # Plot 1: UMAP colored by tumor/normal
  tryCatch({
    p1 <- DimPlot(sobj_labeled, group.by = "class", reduction = "umap",
                  cols = c("tumor" = "#E41A1C", "normal" = "#4DAF4A"),
                  pt.size = 0.2) +
      labs(title = "SCEVAN: Tumor vs Normal") +
      theme_minimal()

    p2 <- DimPlot(sobj_labeled, group.by = "class", reduction = "umap",
                  split.by = "group",
                  cols = c("tumor" = "#E41A1C", "normal" = "#4DAF4A"),
                  pt.size = 0.2) +
      labs(title = "SCEVAN by Group") +
      theme_minimal()

    pdf(file.path(plot_out, "01_scevan_umap_tumor_normal.pdf"), width = 16, height = 6)
    print(p1 | p2)
    dev.off()

    png(file.path(plot_out, "01_scevan_umap_tumor_normal.png"), width = 1600, height = 600, res = 150)
    print(p1 | p2)
    dev.off()

    cat("  ✅ 01_scevan_umap_tumor_normal\n")
  }, error = function(e) cat(sprintf("  ⚠️ UMAP plot failed: %s\n", e$message)))

  # Plot 2: UMAP split by sample, colored by class
  tryCatch({
    p3 <- DimPlot(sobj_labeled, group.by = "class", reduction = "umap",
                  split.by = "sample_id",
                  cols = c("tumor" = "#E41A1C", "normal" = "#4DAF4A"),
                  pt.size = 0.1, ncol = 3) +
      labs(title = "SCEVAN by Sample") +
      theme_minimal()

    pdf(file.path(plot_out, "02_scevan_umap_by_sample.pdf"), width = 18, height = 10)
    print(p3)
    dev.off()

    png(file.path(plot_out, "02_scevan_umap_by_sample.png"), width = 1800, height = 1000, res = 150)
    print(p3)
    dev.off()

    cat("  ✅ 02_scevan_umap_by_sample\n")
  }, error = function(e) cat(sprintf("  ⚠️ Sample UMAP failed: %s\n", e$message)))

  # Plot 3: Tumor fraction bar chart
  tryCatch({
    p4 <- ggplot(cluster_summary, aes(x = reorder(cluster, -pct_tumor), y = pct_tumor)) +
      geom_col(aes(fill = pct_tumor > 50), show.legend = FALSE) +
      scale_fill_manual(values = c("FALSE" = "#4DAF4A", "TRUE" = "#E41A1C")) +
      geom_hline(yintercept = 50, linetype = "dashed", color = "grey40") +
      geom_text(aes(label = sprintf("%.0f%%", pct_tumor)), vjust = -0.3, size = 3) +
      labs(title = "Tumor Fraction by Cluster",
           x = "Cluster", y = "% Tumor") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 0))

    p5 <- ggplot(group_summary, aes(x = group, y = pct_tumor, fill = group)) +
      geom_col(show.legend = FALSE, alpha = 0.8) +
      geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct_tumor, n_tumor)), vjust = -0.1, size = 3.5) +
      labs(title = "Tumor Fraction by Group", x = "", y = "% Tumor") +
      theme_minimal()

    pdf(file.path(plot_out, "03_scevan_tumor_fraction.pdf"), width = 14, height = 6)
    print(p4 | p5)
    dev.off()

    png(file.path(plot_out, "03_scevan_tumor_fraction.png"), width = 1400, height = 600, res = 150)
    print(p4 | p5)
    dev.off()

    cat("  ✅ 03_scevan_tumor_fraction\n")
  }, error = function(e) cat(sprintf("  ⚠️ Bar chart failed: %s\n", e$message)))

} else {
  cat("⚠️  No class column — skipping summaries and plots\n")
  saveRDS(results, file.path(main_out, "scevan_results_raw.rds"))
}

# ============================================================================
# 11. Copy SCEVAN-generated files
# ============================================================================

scevan_files <- list.files(local_out, pattern = "\\.(pdf|png|csv|txt)$",
                           full.names = TRUE, recursive = TRUE)
if (length(scevan_files) > 0) {
  # Copy SCEVAN's own plots to main plot dir
  scevan_plots <- grep("\\.(pdf|png)$", scevan_files, value = TRUE)
  if (length(scevan_plots) > 0) {
    file.copy(scevan_plots, plot_out, overwrite = TRUE)
    cat(sprintf("\n  %d SCEVAN-generated plot(s) → %s\n", length(scevan_plots), plot_out))
  }
}

# ============================================================================
# 12. Final Report
# ============================================================================

gc()

cat("\n")
cat("╔══════════════════════════════════════════════════════╗\n")
cat("║                 SCEVAN Complete                      ║\n")
cat(sprintf("║  Elapsed:     %.1f min                               ║\n", elapsed))
cat(sprintf("║  Cells:       %-6d                                  ║\n", ncol(sobj)))
if (!is.null(class_col)) {
  n_tum <- sum(scevan_df[[class_col]] == "tumor", na.rm = TRUE)
  n_nor <- sum(scevan_df[[class_col]] == "normal", na.rm = TRUE)
  cat(sprintf("║  Tumor:       %-6d (%.1f%%)                         ║\n",
              n_tum, n_tum / nrow(scevan_df) * 100))
  cat(sprintf("║  Normal:      %-6d (%.1f%%)                         ║\n",
              n_nor, n_nor / nrow(scevan_df) * 100))
}
cat("╠══════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Module:  modules/cnv/results/scevan/                ║\n"))
cat(sprintf("║  Main:    results/scrna/10_cnv/                      ║\n"))
cat(sprintf("║  Plots:   results/scrna/10_cnv/plots/               ║\n"))
cat(sprintf("║  Reports: results/scrna/10_cnv/reports/             ║\n"))
cat("╚══════════════════════════════════════════════════════╝\n")
