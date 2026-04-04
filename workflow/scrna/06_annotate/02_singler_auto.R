#!/usr/bin/env Rscript
# ============================================================================
# 02_singler_auto.R — Automatic cell type annotation using SingleR
#
# References:
#   - celldex::MouseRNAseqData()  (broad cell types)
#   - celldex::ImmGenData()       (immune-focused)
#
# Runs both per-cell and per-cluster annotation
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(SingleCellExperiment)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

source(here("workflow", "scrna", "functions", "io_scrna.R"))
source(here("workflow", "scrna", "functions", "qc_utils.R"))
source(here("scripts", "utils", "utils_plotting.R"))

# ---- Config ----
out_base <- here("results", "scrna", "06_annotate")
dirs <- list(
  singler = file.path(out_base, "singler"),
  plots   = file.path(out_base, "plots"),
  reports = file.path(out_base, "reports")
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat("==============================================================\n")
cat("   Phase 6 Step 02: SingleR Auto Annotation                   \n")
cat("==============================================================\n\n")

# ---- Load ----
log_msg("Loading clustered object...")
sobj <- readRDS(here("results", "scrna", "05_cluster", "objects", "seurat_clustered.rds"))
DefaultAssay(sobj) <- "SCT"
Idents(sobj) <- "seurat_clusters"
log_msg(sprintf("  %d cells, %d clusters", ncol(sobj), length(levels(Idents(sobj)))))

# ---- Convert to SCE ----
log_msg("Converting to SCE...")
sce <- as.SingleCellExperiment(sobj)

# ---- Load references ----
log_msg("Loading MouseRNAseqData...")
ref_mouse <- celldex::MouseRNAseqData()
log_msg(sprintf("  %d genes, %d samples, %d labels",
                nrow(ref_mouse), ncol(ref_mouse), length(unique(ref_mouse$label.main))))

log_msg("Loading ImmGenData...")
ref_immgen <- celldex::ImmGenData()
log_msg(sprintf("  %d genes, %d samples, %d labels",
                nrow(ref_immgen), ncol(ref_immgen), length(unique(ref_immgen$label.main))))

# ============================================================================
# Per-cluster annotation (robust, fast, recommended for initial annotation)
# ============================================================================

log_msg("SingleR per-cluster: MouseRNAseqData...")
t0 <- Sys.time()
pred_mouse_cluster <- SingleR(
  test      = sce,
  ref       = ref_mouse,
  labels    = ref_mouse$label.main,
  clusters  = sobj$seurat_clusters,
  de.method = "wilcox"
)
log_msg(sprintf("  Done in %.1f min", difftime(Sys.time(), t0, units = "mins")))

log_msg("SingleR per-cluster: ImmGenData...")
t0 <- Sys.time()
pred_immgen_cluster <- SingleR(
  test      = sce,
  ref       = ref_immgen,
  labels    = ref_immgen$label.main,
  clusters  = sobj$seurat_clusters,
  de.method = "wilcox"
)
log_msg(sprintf("  Done in %.1f min", difftime(Sys.time(), t0, units = "mins")))

# ============================================================================
# Per-cell annotation (slower, but useful for checking heterogeneity)
# ============================================================================

log_msg("SingleR per-cell: MouseRNAseqData...")
t0 <- Sys.time()
pred_mouse_cell <- SingleR(
  test      = sce,
  ref       = ref_mouse,
  labels    = ref_mouse$label.main,
  de.method = "wilcox"
)
log_msg(sprintf("  Done in %.1f min", difftime(Sys.time(), t0, units = "mins")))

log_msg("SingleR per-cell: ImmGenData...")
t0 <- Sys.time()
pred_immgen_cell <- SingleR(
  test      = sce,
  ref       = ref_immgen,
  labels    = ref_immgen$label.main,
  de.method = "wilcox"
)
log_msg(sprintf("  Done in %.1f min", difftime(Sys.time(), t0, units = "mins")))

# ---- Save SingleR objects ----
saveRDS(pred_mouse_cluster,  file.path(dirs$singler, "singler_MouseRNAseq_percluster.rds"))
saveRDS(pred_immgen_cluster, file.path(dirs$singler, "singler_ImmGen_percluster.rds"))
saveRDS(pred_mouse_cell,     file.path(dirs$singler, "singler_MouseRNAseq_percell.rds"))
saveRDS(pred_immgen_cell,    file.path(dirs$singler, "singler_ImmGen_percell.rds"))

# ---- Add to Seurat metadata ----
# Per-cell labels
#sobj$singler_mouse_cell   <- pred_mouse_cell$labels
#sobj$singler_immgen_cell  <- pred_immgen_cell$labels

# 新（绕过 Seurat name matching）：
sobj@meta.data$singler_mouse_cell  <- pred_mouse_cell$labels[match(colnames(sobj), rownames(pred_mouse_cell))]
sobj@meta.data$singler_immgen_cell <- pred_immgen_cell$labels[match(colnames(sobj), rownames(pred_immgen_cell))]


# Per-cluster labels (map back to cells)
cl_labels_mouse  <- pred_mouse_cluster$labels
cl_labels_immgen <- pred_immgen_cluster$labels
names(cl_labels_mouse)  <- rownames(pred_mouse_cluster)
names(cl_labels_immgen) <- rownames(pred_immgen_cluster)

sobj@meta.data$singler_mouse_cluster  <- cl_labels_mouse[as.character(sobj$seurat_clusters)]
sobj@meta.data$singler_immgen_cluster <- cl_labels_immgen[as.character(sobj$seurat_clusters)]

# ---- Consensus table ----
log_msg("Building consensus annotation table...")

cluster_anno <- data.frame(
  cluster          = rownames(pred_mouse_cluster),
  mouse_label      = pred_mouse_cluster$labels,
  mouse_score      = round(apply(pred_mouse_cluster$scores, 1, max), 3),
  mouse_pruned     = pred_mouse_cluster$pruned.labels,
  immgen_label     = pred_immgen_cluster$labels,
  immgen_score     = round(apply(pred_immgen_cluster$scores, 1, max), 3),
  immgen_pruned    = pred_immgen_cluster$pruned.labels,
  stringsAsFactors = FALSE
)

# Add cell counts
cluster_sizes <- table(sobj$seurat_clusters)
cluster_anno$n_cells <- as.integer(cluster_sizes[cluster_anno$cluster])
cluster_anno$refs_agree <- cluster_anno$mouse_label == cluster_anno$immgen_label
cluster_anno <- cluster_anno %>% arrange(as.integer(cluster))

fwrite(cluster_anno, file.path(dirs$reports, "singler_cluster_annotation.csv"))
cat("\n=== SingleR Per-Cluster Annotation ===\n")
print(as.data.frame(cluster_anno), row.names = FALSE)

# ---- Per-cell label distribution within clusters ----
cell_label_dist <- sobj@meta.data %>%
  group_by(seurat_clusters, singler_mouse_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  arrange(seurat_clusters, desc(pct))

fwrite(cell_label_dist, file.path(dirs$reports, "percell_label_distribution.csv"))

# ============================================================================
# Plots
# ============================================================================

# ---- UMAP: cluster vs SingleR labels ----
log_msg("Plotting SingleR UMAPs...")
tryCatch({
  p_cluster <- DimPlot(sobj, group.by = "seurat_clusters", label = TRUE, repel = TRUE,
                       pt.size = 0.1, label.size = 3) +
    labs(title = "Clusters (res=0.8)") + NoLegend()

  p_mouse <- DimPlot(sobj, group.by = "singler_mouse_cluster", label = TRUE, repel = TRUE,
                     pt.size = 0.1, label.size = 3) +
    labs(title = "SingleR: MouseRNAseqData") +
    theme(legend.text = element_text(size = 7))

  p_immgen <- DimPlot(sobj, group.by = "singler_immgen_cluster", label = TRUE, repel = TRUE,
                      pt.size = 0.1, label.size = 3) +
    labs(title = "SingleR: ImmGenData") +
    theme(legend.text = element_text(size = 7))

  pdf(file.path(dirs$plots, "06_singler_umap.pdf"), width = 20, height = 16)
  print((p_cluster | p_mouse) / (p_cluster | p_immgen))
  dev.off()
  png(file.path(dirs$plots, "06_singler_umap.png"), width = 2000, height = 1600, res = 150)
  print((p_cluster | p_mouse) / (p_cluster | p_immgen))
  dev.off()
  log_msg("  Done: 06_singler_umap")
}, error = function(e) log_msg(sprintf("  UMAP failed: %s", e$message), "warn"))

# ---- Score heatmaps ----
log_msg("Plotting score heatmaps...")
tryCatch({
  pdf(file.path(dirs$plots, "07_singler_scores_mouse.pdf"), width = 12, height = 8)
  plotScoreHeatmap(pred_mouse_cluster, show.pruned = TRUE,
                   main = "SingleR Scores: MouseRNAseqData (per-cluster)")
  dev.off()
  pdf(file.path(dirs$plots, "07_singler_scores_immgen.pdf"), width = 14, height = 8)
  plotScoreHeatmap(pred_immgen_cluster, show.pruned = TRUE,
                   main = "SingleR Scores: ImmGenData (per-cluster)")
  dev.off()
  log_msg("  Done: 07_singler_scores")
}, error = function(e) log_msg(sprintf("  Score heatmap failed: %s", e$message), "warn"))

# ---- Delta distribution ----
tryCatch({
  pdf(file.path(dirs$plots, "08_singler_delta_mouse.pdf"), width = 14, height = 6)
  plotDeltaDistribution(pred_mouse_cluster,
                        main = "Delta Distribution: MouseRNAseqData")
  dev.off()
  pdf(file.path(dirs$plots, "08_singler_delta_immgen.pdf"), width = 14, height = 6)
  plotDeltaDistribution(pred_immgen_cluster,
                        main = "Delta Distribution: ImmGenData")
  dev.off()
  log_msg("  Done: 08_singler_delta")
}, error = function(e) log_msg(sprintf("  Delta plot failed: %s", e$message), "warn"))

# ---- Save annotated object ----
log_msg("Saving annotated object...")
saveRDS(sobj, file.path(out_base, "objects", "seurat_singler_annotated.rds"))
log_msg(sprintf("  -> %s", file.path(out_base, "objects", "seurat_singler_annotated.rds")))

rm(sce); gc()
cat("\n========== SingleR annotation complete ==========\n")
