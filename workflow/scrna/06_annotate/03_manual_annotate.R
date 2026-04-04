#!/usr/bin/env Rscript
# ============================================================================
# 03_manual_annotate.R — Apply manual annotation from celltype_mapping.csv
#
# Input:
#   1. seurat_clustered.rds          (Phase 5)
#   2. configs/annotation/celltype_mapping.csv  (human-curated)
#   3. SingleR results               (Phase 6 Step 2, optional)
#
# Output:
#   1. seurat_annotated_full.rds     (all cells, artifact labeled but kept)
#   2. seurat_annotated.rds          (artifact removed, re-embedded)
#   3. Plots: annotated UMAPs, composition bars, QC summary
#   4. Reports: cell counts per annotation, per group
#
# Design: This script contains ZERO hardcoded cell type names.
#         All annotation comes from the CSV.
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(data.table)
  library(scales)
})

source(here("scripts", "utils", "utils_io.R"))
source(here("scripts", "utils", "utils_plotting.R"))
source(here("workflow", "scrna", "functions", "qc_utils.R"))

# ---- Paths ----
input_obj   <- here("results", "scrna", "05_cluster", "objects", "seurat_clustered.rds")
mapping_csv <- here("configs", "annotation", "celltype_mapping.csv")
out_base    <- here("results", "scrna", "06_annotate")
obj_dir     <- file.path(out_base, "objects")
plot_dir    <- file.path(out_base, "plots")
report_dir  <- file.path(out_base, "reports")
for (d in c(obj_dir, plot_dir, report_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

cat("\n")
cat("╔══════════════════════════════════════════════════╗\n")
cat("║   Phase 6 Step 03: Manual Annotation             ║\n")
cat("╚══════════════════════════════════════════════════╝\n\n")

# ============================================================================
# 1. Load inputs
# ============================================================================
log_msg("Loading clustered object...")
sobj <- readRDS(input_obj)
log_msg(sprintf("  %d cells, %d clusters", ncol(sobj), length(unique(sobj$seurat_clusters))))

log_msg("Loading celltype_mapping.csv...")
if (!file.exists(mapping_csv)) stop("celltype_mapping.csv not found: ", mapping_csv)
mapping <- fread(mapping_csv)
log_msg(sprintf("  %d rows in mapping", nrow(mapping)))

# Validate mapping
required_cols <- c("cluster", "celltype_L1", "celltype_L2", "confidence", "key_markers", "notes")
missing_cols <- setdiff(required_cols, colnames(mapping))
if (length(missing_cols) > 0) stop("Missing columns in mapping CSV: ", paste(missing_cols, collapse = ", "))

# Check all clusters are mapped
clusters_in_data <- sort(as.integer(as.character(unique(sobj$seurat_clusters))))
clusters_in_csv  <- sort(mapping$cluster)
unmapped <- setdiff(clusters_in_data, clusters_in_csv)
if (length(unmapped) > 0) {
  log_msg(sprintf("⚠️  Unmapped clusters: %s", paste(unmapped, collapse = ", ")), "warn")
}

# ============================================================================
# 2. Write annotations to metadata
# ============================================================================
log_msg("Writing annotations to metadata...")

# Create lookup vectors
cl_to_L1 <- setNames(mapping$celltype_L1, as.character(mapping$cluster))
cl_to_L2 <- setNames(mapping$celltype_L2, as.character(mapping$cluster))
cl_to_conf <- setNames(mapping$confidence, as.character(mapping$cluster))

cl_char <- as.character(sobj$seurat_clusters)

sobj@meta.data$celltype_L1  <- unname(cl_to_L1[cl_char])
sobj@meta.data$celltype_L2  <- unname(cl_to_L2[cl_char])
sobj@meta.data$annotation_confidence <- unname(cl_to_conf[cl_char])

# Mark artifacts
artifact_clusters <- mapping$cluster[mapping$celltype_L1 == "Artifact"]
sobj@meta.data$is_artifact <- cl_char %in% as.character(artifact_clusters)

n_artifact <- sum(sobj$is_artifact)
log_msg(sprintf("  Artifact cells: %d (%.1f%%) — clusters: %s",
                n_artifact, n_artifact / ncol(sobj) * 100,
                paste(artifact_clusters, collapse = ", ")))

# ============================================================================
# 3. Save FULL annotated object (artifacts labeled but NOT removed)
# ============================================================================
log_msg("Saving full annotated object (with artifacts)...")
saveRDS(sobj, file.path(obj_dir, "seurat_annotated_full.rds"))

# ============================================================================
# 4. Remove artifacts + re-embed UMAP
# ============================================================================
log_msg("Removing artifact clusters and re-embedding...")
sobj_clean <- subset(sobj, is_artifact == FALSE)
log_msg(sprintf("  %d → %d cells after artifact removal", ncol(sobj), ncol(sobj_clean)))

# Re-run UMAP on clean cells (using existing Harmony reduction)
log_msg("  Re-running UMAP on clean cells...")
sobj_clean <- RunUMAP(sobj_clean, reduction = "harmony", dims = 1:30,
                       verbose = FALSE)

# Set default identity to L2
Idents(sobj_clean) <- "celltype_L2"

# Save clean annotated object
saveRDS(sobj_clean, file.path(obj_dir, "seurat_annotated.rds"))
log_msg(sprintf("  ✅ seurat_annotated.rds: %d cells", ncol(sobj_clean)))

# ============================================================================
# 5. Reports
# ============================================================================
log_msg("Generating reports...")

# 5a. Cell count per L1/L2
count_L1 <- sobj_clean@meta.data %>%
  group_by(celltype_L1) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  mutate(pct = round(n_cells / sum(n_cells) * 100, 2)) %>%
  arrange(desc(n_cells))
fwrite(count_L1, file.path(report_dir, "cellcount_by_L1.csv"))

count_L2 <- sobj_clean@meta.data %>%
  group_by(celltype_L1, celltype_L2, seurat_clusters) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  mutate(pct = round(n_cells / sum(n_cells) * 100, 2)) %>%
  arrange(celltype_L1, desc(n_cells))
fwrite(count_L2, file.path(report_dir, "cellcount_by_L2.csv"))

# 5b. Cell count per L1 × group
count_L1_group <- sobj_clean@meta.data %>%
  group_by(group, celltype_L1) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(pct = round(n_cells / sum(n_cells) * 100, 2)) %>%
  arrange(group, desc(n_cells))
fwrite(count_L1_group, file.path(report_dir, "cellcount_L1_by_group.csv"))

count_L2_group <- sobj_clean@meta.data %>%
  group_by(group, celltype_L1, celltype_L2) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(pct = round(n_cells / sum(n_cells) * 100, 2)) %>%
  arrange(group, desc(n_cells))
fwrite(count_L2_group, file.path(report_dir, "cellcount_L2_by_group.csv"))

# 5c. Annotation summary
anno_summary <- mapping %>%
  left_join(
    sobj_clean@meta.data %>%
      group_by(seurat_clusters) %>%
      summarise(n_cells_clean = n(), .groups = "drop") %>%
      mutate(seurat_clusters = as.integer(as.character(seurat_clusters))),
    by = c("cluster" = "seurat_clusters")
  )
anno_summary$n_cells_clean[is.na(anno_summary$n_cells_clean)] <- 0
fwrite(anno_summary, file.path(report_dir, "annotation_summary.csv"))

cat("\n=== Cell Counts by L1 ===\n")
print(as.data.frame(count_L1), row.names = FALSE)

# ============================================================================
# 6. Plots
# ============================================================================
log_msg("Generating annotation plots...")

# ---- Color palette ----
# Deterministic palette: one color per L1, shades per L2
l1_levels <- sort(unique(sobj_clean$celltype_L1))
n_l1 <- length(l1_levels)

# Base colors for L1 categories
l1_base_colors <- c(
  "Cycling"          = "#FFD700",
  "Epithelial_normal"= "#2ca02c",
  "Lymphoid"         = "#1f77b4",
  "Myeloid"          = "#ff7f0e",
  "Stromal"          = "#9467bd",
  "Tumor_putative"   = "#d62728",
  "Unresolved"       = "#7f7f7f"
)
# Fallback for any unexpected L1
for (l in l1_levels) {
  if (is.na(l1_base_colors[l])) l1_base_colors[l] <- "grey50"
}

# L2 palette: generate shades within each L1
l2_colors <- c()
for (l1 in l1_levels) {
  l2s <- sort(unique(sobj_clean$celltype_L2[sobj_clean$celltype_L1 == l1]))
  n <- length(l2s)
  if (n == 1) {
    cols <- l1_base_colors[l1]
  } else {
    base <- col2rgb(l1_base_colors[l1]) / 255
    cols <- colorRampPalette(c(
      rgb(min(base[1]+0.3,1), min(base[2]+0.3,1), min(base[3]+0.3,1)),
      l1_base_colors[l1],
      rgb(max(base[1]-0.2,0), max(base[2]-0.2,0), max(base[3]-0.2,0))
    ))(n)
  }
  names(cols) <- l2s
  l2_colors <- c(l2_colors, cols)
}

# ---- Plot 10: UMAP by L1 ----
p_l1 <- DimPlot(sobj_clean, group.by = "celltype_L1", label = TRUE, repel = TRUE,
                pt.size = 0.2, label.size = 3.5) +
  scale_color_manual(values = l1_base_colors) +
  labs(title = "Cell Type (L1)") +
  theme_project() + NoLegend()

p_l1_leg <- DimPlot(sobj_clean, group.by = "celltype_L1", label = FALSE,
                    pt.size = 0.2) +
  scale_color_manual(values = l1_base_colors) +
  labs(title = "Cell Type (L1)") +
  theme_project()

pdf(file.path(plot_dir, "10_umap_celltype_L1.pdf"), width = 16, height = 7)
print(p_l1 | p_l1_leg)
dev.off()
png(file.path(plot_dir, "10_umap_celltype_L1.png"), width = 1600, height = 700, res = 150)
print(p_l1 | p_l1_leg)
dev.off()
log_msg("  ✅ 10_umap_celltype_L1")

# ---- Plot 11: UMAP by L2 ----
p_l2 <- DimPlot(sobj_clean, group.by = "celltype_L2", label = TRUE, repel = TRUE,
                pt.size = 0.2, label.size = 2.5) +
  scale_color_manual(values = l2_colors) +
  labs(title = "Cell Type (L2)") +
  theme_project() + NoLegend()

p_l2_leg <- DimPlot(sobj_clean, group.by = "celltype_L2", label = FALSE,
                    pt.size = 0.2) +
  scale_color_manual(values = l2_colors) +
  labs(title = "Cell Type (L2)") +
  theme_project() +
  theme(legend.text = element_text(size = 6))

pdf(file.path(plot_dir, "11_umap_celltype_L2.pdf"), width = 20, height = 8)
print(p_l2 | p_l2_leg)
dev.off()
png(file.path(plot_dir, "11_umap_celltype_L2.png"), width = 2000, height = 800, res = 150)
print(p_l2 | p_l2_leg)
dev.off()
log_msg("  ✅ 11_umap_celltype_L2")

# ---- Plot 12: UMAP by group (split) ----
p_grp <- DimPlot(sobj_clean, group.by = "celltype_L1", split.by = "group",
                 label = TRUE, repel = TRUE, pt.size = 0.15, label.size = 2.5) +
  scale_color_manual(values = l1_base_colors) +
  labs(title = "Cell Type L1 by Group") +
  theme_project() + NoLegend()

pdf(file.path(plot_dir, "12_umap_L1_split_by_group.pdf"), width = 20, height = 7)
print(p_grp)
dev.off()
png(file.path(plot_dir, "12_umap_L1_split_by_group.png"), width = 2000, height = 700, res = 150)
print(p_grp)
dev.off()
log_msg("  ✅ 12_umap_L1_split_by_group")

# ---- Plot 13: UMAP original clusters (for cross-reference) ----
p_cl <- DimPlot(sobj_clean, group.by = "seurat_clusters", label = TRUE, repel = TRUE,
                pt.size = 0.2, label.size = 3.5) +
  labs(title = "Original Clusters (post-artifact removal)") +
  theme_project() + NoLegend()

pdf(file.path(plot_dir, "13_umap_clusters_clean.pdf"), width = 10, height = 8)
print(p_cl)
dev.off()
png(file.path(plot_dir, "13_umap_clusters_clean.png"), width = 1000, height = 800, res = 150)
print(p_cl)
dev.off()
log_msg("  ✅ 13_umap_clusters_clean")

# ---- Plot 14: Cell composition stacked bars ----
comp_data <- sobj_clean@meta.data %>%
  group_by(group, celltype_L1) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(pct = n / sum(n) * 100)

p_comp <- ggplot(comp_data, aes(x = group, y = pct, fill = celltype_L1)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(pct > 3, sprintf("%.0f%%", pct), "")),
            position = position_stack(vjust = 0.5), size = 2.5) +
  scale_fill_manual(values = l1_base_colors) +
  labs(title = "Cell Composition by Group (L1)", x = "", y = "Percentage", fill = "Cell Type") +
  theme_project() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(plot_dir, "14_composition_L1_by_group.pdf"), width = 10, height = 8)
print(p_comp)
dev.off()
png(file.path(plot_dir, "14_composition_L1_by_group.png"), width = 1000, height = 800, res = 150)
print(p_comp)
dev.off()
log_msg("  ✅ 14_composition_L1_by_group")

# ---- Plot 15: Cell composition L2 stacked bars ----
comp_l2 <- sobj_clean@meta.data %>%
  group_by(group, celltype_L2) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(pct = n / sum(n) * 100)

p_comp_l2 <- ggplot(comp_l2, aes(x = group, y = pct, fill = celltype_L2)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = l2_colors) +
  labs(title = "Cell Composition by Group (L2)", x = "", y = "Percentage", fill = "Cell Type") +
  theme_project() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"))

pdf(file.path(plot_dir, "15_composition_L2_by_group.pdf"), width = 14, height = 8)
print(p_comp_l2)
dev.off()
png(file.path(plot_dir, "15_composition_L2_by_group.png"), width = 1400, height = 800, res = 150)
print(p_comp_l2)
dev.off()
log_msg("  ✅ 15_composition_L2_by_group")

# ---- Plot 16: Cell count per sample (stacked by L1) ----
comp_sample <- sobj_clean@meta.data %>%
  group_by(sample_id, celltype_L1) %>%
  summarise(n = n(), .groups = "drop")

p_sample <- ggplot(comp_sample, aes(x = sample_id, y = n, fill = celltype_L1)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = l1_base_colors) +
  scale_y_continuous(labels = comma) +
  labs(title = "Cell Count per Sample (L1)", x = "", y = "Cells", fill = "Cell Type") +
  theme_project() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(plot_dir, "16_cellcount_per_sample_L1.pdf"), width = 12, height = 7)
print(p_sample)
dev.off()
png(file.path(plot_dir, "16_cellcount_per_sample_L1.png"), width = 1200, height = 700, res = 150)
print(p_sample)
dev.off()
log_msg("  ✅ 16_cellcount_per_sample_L1")

# ---- Plot 17: Confidence overlay on UMAP ----
conf_colors <- c("high" = "#2ca02c", "medium" = "#FFD700", "low" = "#ff7f0e", "extreme" = "#d62728")
p_conf <- DimPlot(sobj_clean, group.by = "annotation_confidence", pt.size = 0.2) +
  scale_color_manual(values = conf_colors) +
  labs(title = "Annotation Confidence") +
  theme_project()

pdf(file.path(plot_dir, "17_umap_confidence.pdf"), width = 10, height = 8)
print(p_conf)
dev.off()
png(file.path(plot_dir, "17_umap_confidence.png"), width = 1000, height = 800, res = 150)
print(p_conf)
dev.off()
log_msg("  ✅ 17_umap_confidence")

# ============================================================================
# Done
# ============================================================================
gc()

cat("\n")
cat("╔══════════════════════════════════════════════════╗\n")
cat("║   Manual Annotation Complete                     ║\n")
cat("╠══════════════════════════════════════════════════╣\n")
cat(sprintf("║  Total cells (full):   %6d                     ║\n", ncol(sobj)))
cat(sprintf("║  Artifact removed:     %6d                     ║\n", n_artifact))
cat(sprintf("║  Clean cells:          %6d                     ║\n", ncol(sobj_clean)))
cat(sprintf("║  Unique L1 types:      %6d                     ║\n", length(unique(sobj_clean$celltype_L1))))
cat(sprintf("║  Unique L2 types:      %6d                     ║\n", length(unique(sobj_clean$celltype_L2))))
cat("╠══════════════════════════════════════════════════╣\n")
cat(sprintf("║  Objects → %s\n", obj_dir))
cat(sprintf("║  Plots   → %s\n", plot_dir))
cat(sprintf("║  Reports → %s\n", report_dir))
cat("╚══════════════════════════════════════════════════╝\n\n")
