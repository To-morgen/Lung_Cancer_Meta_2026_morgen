#!/usr/bin/env Rscript
# ============================================================================
# 01_find_markers.R — FindAllMarkers + lineage analysis + project gene analysis
#
# ALL gene lists and parameters come from:
#   configs/params/scrna_annotation_params.yaml
#
# This script contains ZERO hardcoded gene names.
# To change markers/genes: edit the YAML, not this file.
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

source(here("workflow", "scrna", "functions", "io_scrna.R"))
source(here("workflow", "scrna", "functions", "qc_utils.R"))
source(here("workflow", "scrna", "functions", "annotation_utils.R"))
source(here("scripts", "utils", "utils_plotting.R"))

# ---- Config (from YAML, not hardcoded) ----
ANNO_PARAMS <- load_annotation_params()
out_base    <- here("results", "scrna", "06_annotate")
dirs <- list(
  markers = file.path(out_base, "markers"),
  plots   = file.path(out_base, "plots"),
  reports = file.path(out_base, "reports")
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat("==============================================================\n")
cat("   Phase 6 Step 01: FindAllMarkers                            \n")
cat("==============================================================\n\n")

# ---- Load ----
log_msg("Loading clustered object...")
sobj <- readRDS(here("results", "scrna", "05_cluster", "objects", "seurat_clustered.rds"))
DefaultAssay(sobj) <- "SCT"
Idents(sobj) <- "seurat_clusters"
n_clusters <- length(levels(Idents(sobj)))
log_msg(sprintf("  %d cells, %d clusters", ncol(sobj), n_clusters))

# ---- PrepSCTFindMarkers ----
log_msg("Running PrepSCTFindMarkers...")
sobj <- PrepSCTFindMarkers(sobj)

# ---- FindAllMarkers (params from YAML) ----
fm_params <- ANNO_PARAMS$find_markers
log_msg(sprintf("FindAllMarkers: test=%s, min_pct=%.2f, logfc=%.2f, only_pos=%s",
                fm_params$test_use, fm_params$min_pct,
                fm_params$logfc_threshold, fm_params$only_pos))

t0 <- Sys.time()
markers <- FindAllMarkers(
  sobj,
  only.pos        = fm_params$only_pos,
  min.pct         = fm_params$min_pct,
  logfc.threshold = fm_params$logfc_threshold,
  test.use        = fm_params$test_use,
  verbose         = TRUE
)
elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
log_msg(sprintf("  Done in %.1f min — %d markers", elapsed, nrow(markers)))

# ---- Save full + top-N tables ----
fwrite(markers, file.path(dirs$markers, "all_markers.csv"))

for (n in fm_params$top_n) {
  top <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = n)
  fwrite(top, file.path(dirs$markers, sprintf("top%d_markers_per_cluster.csv", n)))
}

# ---- Summary ----
marker_count <- markers %>%
  group_by(cluster) %>%
  summarise(
    n_markers    = n(),
    n_sig        = sum(p_val_adj < 0.05),
    top_gene     = gene[which.max(avg_log2FC)],
    top_logfc    = round(max(avg_log2FC), 2),
    median_logfc = round(median(avg_log2FC), 2),
    .groups      = "drop"
  ) %>%
  arrange(as.integer(as.character(cluster)))

fwrite(marker_count, file.path(dirs$reports, "marker_summary_per_cluster.csv"))
cat("\n=== Markers per Cluster ===\n")
print(as.data.frame(marker_count), row.names = FALSE)

# ============================================================================
# Plot 01: Heatmap (top 5)
# ============================================================================
log_msg("Plotting marker heatmap...")
set.seed(42)
sobj_ds <- subset(sobj, downsample = 200)
top5 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)

tryCatch({
  pdf(file.path(dirs$plots, "01_marker_heatmap_top5.pdf"), width = 18, height = 14)
  print(DoHeatmap(sobj_ds, features = top5$gene, size = 2.5, angle = 90) +
          theme(axis.text.y = element_text(size = 5)) +
          labs(title = "Top 5 Markers per Cluster"))
  dev.off()
  png(file.path(dirs$plots, "01_marker_heatmap_top5.png"), width = 1800, height = 1400, res = 150)
  print(DoHeatmap(sobj_ds, features = top5$gene, size = 2.5, angle = 90) +
          theme(axis.text.y = element_text(size = 5)) +
          labs(title = "Top 5 Markers per Cluster"))
  dev.off()
  log_msg("  Done: 01_marker_heatmap_top5")
}, error = function(e) log_msg(sprintf("  Heatmap failed: %s", e$message), "warn"))
rm(sobj_ds); gc()

# ============================================================================
# Plot 02: Dotplot (top 3)
# ============================================================================
log_msg("Plotting top3 dotplot...")
top3 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 3)

tryCatch({
  pdf(file.path(dirs$plots, "02_marker_dotplot_top3.pdf"), width = 20, height = 10)
  print(DotPlot(sobj, features = unique(top3$gene), cluster.idents = TRUE) +
          coord_flip() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                axis.text.y = element_text(size = 6)) +
          labs(title = "Top 3 Markers per Cluster"))
  dev.off()
  png(file.path(dirs$plots, "02_marker_dotplot_top3.png"), width = 2000, height = 1000, res = 150)
  print(DotPlot(sobj, features = unique(top3$gene), cluster.idents = TRUE) +
          coord_flip() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                axis.text.y = element_text(size = 6)) +
          labs(title = "Top 3 Markers per Cluster"))
  dev.off()
  log_msg("  Done: 02_marker_dotplot_top3")
}, error = function(e) log_msg(sprintf("  Dotplot failed: %s", e$message), "warn"))

# ============================================================================
# Plot 03: Lineage dotplot (from YAML config, deduplicated)
# ============================================================================
log_msg("Plotting lineage markers (from config)...")

lineage_genes <- get_lineage_markers(ANNO_PARAMS)  # auto-deduplicated
lineage_ref   <- get_lineage_reference(ANNO_PARAMS)
lineage_ref$present <- lineage_ref$gene %in% rownames(sobj)
fwrite(lineage_ref, file.path(dirs$reports, "lineage_markers_reference.csv"))

lineage_check <- filter_genes_present(lineage_genes, sobj)
log_msg(sprintf("  Lineage: %d present, %d missing", 
                length(lineage_check$present), length(lineage_check$missing)))

if (length(lineage_check$present) > 0) {
  tryCatch({
    pdf(file.path(dirs$plots, "03_lineage_dotplot.pdf"), width = 18, height = 16)
    print(DotPlot(sobj, features = lineage_check$present, cluster.idents = TRUE) +
            coord_flip() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                  axis.text.y = element_text(size = 7)) +
            labs(title = "Lineage Markers (from config, deduplicated)"))
    dev.off()
    png(file.path(dirs$plots, "03_lineage_dotplot.png"), width = 1800, height = 1600, res = 150)
    print(DotPlot(sobj, features = lineage_check$present, cluster.idents = TRUE) +
            coord_flip() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                  axis.text.y = element_text(size = 7)) +
            labs(title = "Lineage Markers (from config, deduplicated)"))
    dev.off()
    log_msg("  Done: 03_lineage_dotplot")
  }, error = function(e) log_msg(sprintf("  Lineage dotplot failed: %s", e$message), "warn"))
}

# ============================================================================
# Plot 04: Key marker FeaturePlots (from YAML config)
# ============================================================================
log_msg("Plotting key feature UMAPs (from config)...")

key_check <- filter_genes_present(ANNO_PARAMS$key_feature_markers, sobj)

if (length(key_check$present) > 0) {
  tryCatch({
    pdf(file.path(dirs$plots, "04_lineage_featureplot.pdf"), width = 20, height = 20)
    print(FeaturePlot(sobj, features = key_check$present, ncol = 4, order = TRUE,
                      cols = c("lightgrey", "darkred"), pt.size = 0.1) &
            theme(plot.title = element_text(size = 10)))
    dev.off()
    png(file.path(dirs$plots, "04_lineage_featureplot.png"), width = 2000, height = 2000, res = 150)
    print(FeaturePlot(sobj, features = key_check$present, ncol = 4, order = TRUE,
                      cols = c("lightgrey", "darkred"), pt.size = 0.1) &
            theme(plot.title = element_text(size = 10)))
    dev.off()
    log_msg("  Done: 04_lineage_featureplot")
  }, error = function(e) log_msg(sprintf("  FeaturePlot failed: %s", e$message), "warn"))
}

# ============================================================================
# Plot 05: Project-specific genes (from YAML config)
# ============================================================================
project_genes <- ANNO_PARAMS$project_genes

if (!is.null(project_genes) && length(project_genes) > 0) {
  log_msg(sprintf("Analyzing %d project-specific gene(s): %s",
                  length(project_genes), paste(project_genes, collapse = ", ")))
  analyze_project_genes(sobj, project_genes, dirs$plots, dirs$reports,
                        prefix = "05_project_gene")
} else {
  log_msg("No project-specific genes defined in config — skipping")
}

gc()
cat("\n========== FindAllMarkers complete ==========\n")
cat(sprintf("Markers -> %s\n", dirs$markers))
cat(sprintf("Plots   -> %s\n", dirs$plots))
cat(sprintf("Reports -> %s\n", dirs$reports))
