#!/usr/bin/env Rscript
# ============================================================================
# 02_singler_auto.R — Automatic cell type annotation using SingleR
#
# Config-driven: references come from annotation params YAML
#   mouse → MouseRNAseqData, ImmGenData
#   human → HumanPrimaryCellAtlasData, BlueprintEncodeData
#
# Runs both per-cluster and (optionally) per-cell annotation
#
# Design: This script contains ZERO hardcoded reference names.
#         All references come from the YAML config via load_annotation_params().
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
  library(qs)
})

source(here("workflow", "scrna", "functions", "io_scrna.R"))
source(here("workflow", "scrna", "functions", "qc_utils.R"))
source(here("workflow", "scrna", "functions", "annotation_utils.R"))
source(here("scripts", "utils", "utils_plotting.R"))

# ---- Config (from YAML, not hardcoded) ----
ANNO_PARAMS <- load_annotation_params()
singler_cfg <- ANNO_PARAMS$singler
run_per_cell <- isTRUE(singler_cfg$run_per_cell)
de_method <- singler_cfg$de_method %||% "wilcox"

out_base <- scrna_base("06_annotate")
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

# ============================================================================
# 1. Load clustered object (qs preferred, rds fallback)
# ============================================================================
log_msg("Loading clustered object...")
input_qs  <- scrna_base("05_cluster", "objects", "seurat_clustered.qs")
input_rds <- scrna_base("05_cluster", "objects", "seurat_clustered.rds")

if (file.exists(input_qs)) {
  sobj <- qread(input_qs)
  log_msg(sprintf("  Loaded from QS: %s", basename(input_qs)))
} else if (file.exists(input_rds)) {
  sobj <- readRDS(input_rds)
  log_msg(sprintf("  Loaded from RDS: %s", basename(input_rds)))
} else {
  stop("No clustered object found in 05_cluster/objects/")
}

DefaultAssay(sobj) <- "SCT"
Idents(sobj) <- "seurat_clusters"
log_msg(sprintf("  %d cells, %d clusters", ncol(sobj), length(levels(Idents(sobj)))))

# ============================================================================
# 2. Convert to SCE
# ============================================================================
log_msg("Converting to SCE...")
sce <- as.SingleCellExperiment(sobj)

# ============================================================================
# 3. Load references dynamically from config
# ============================================================================
refs <- list()
for (ref_cfg in singler_cfg$references) {
  ref_name  <- ref_cfg$name
  ref_pkg   <- ref_cfg$package %||% "celldex"
  label_col <- ref_cfg$label_col %||% "label.main"

  log_msg(sprintf("Loading reference: %s::%s() [labels: %s]", ref_pkg, ref_name, label_col))

  # Dynamically call celldex::<FunctionName>()
  ref_fn <- tryCatch(
    getExportedValue(ref_pkg, ref_name),
    error = function(e) {
      stop(sprintf("Cannot find %s::%s() — is the package installed? Error: %s",
                   ref_pkg, ref_name, e$message))
    }
  )
  ref_data <- ref_fn()
  labels <- ref_data[[label_col]]

  log_msg(sprintf("  %d genes, %d samples, %d unique labels",
                  nrow(ref_data), ncol(ref_data), length(unique(labels))))

  refs[[ref_name]] <- list(
    data      = ref_data,
    labels    = labels,
    label_col = label_col,
    short     = gsub("Data$", "", ref_name)  # e.g. "MouseRNAseq", "ImmGen", "HumanPrimaryCellAtlas"
  )
}

log_msg(sprintf("Loaded %d reference(s): %s", length(refs), paste(names(refs), collapse = ", ")))

# ============================================================================
# 4. Run SingleR per-cluster + (optionally) per-cell for each reference
# ============================================================================
pred_cluster <- list()
pred_cell    <- list()

for (ref_name in names(refs)) {
  ref <- refs[[ref_name]]
  short <- ref$short

  # ── Per-cluster ──
  log_msg(sprintf("SingleR per-cluster: %s...", ref_name))
  t0 <- Sys.time()
  pred_cluster[[ref_name]] <- SingleR(
    test      = sce,
    ref       = ref$data,
    labels    = ref$labels,
    clusters  = sobj$seurat_clusters,
    de.method = de_method
  )
  log_msg(sprintf("  Done in %.1f min", difftime(Sys.time(), t0, units = "mins")))

  # Save
  saveRDS(pred_cluster[[ref_name]],
          file.path(dirs$singler, sprintf("singler_%s_percluster.rds", short)))

  # Map cluster labels back to cells
  cl_labels <- pred_cluster[[ref_name]]$labels
  names(cl_labels) <- rownames(pred_cluster[[ref_name]])
  col_name_cluster <- sprintf("singler_%s_cluster", short)
  sobj@meta.data[[col_name_cluster]] <- cl_labels[as.character(sobj$seurat_clusters)]

  # ── Per-cell (optional, slow) ──
  if (run_per_cell) {
    log_msg(sprintf("SingleR per-cell: %s...", ref_name))
    t0 <- Sys.time()
    pred_cell[[ref_name]] <- SingleR(
      test      = sce,
      ref       = ref$data,
      labels    = ref$labels,
      de.method = de_method
    )
    log_msg(sprintf("  Done in %.1f min", difftime(Sys.time(), t0, units = "mins")))

    saveRDS(pred_cell[[ref_name]],
            file.path(dirs$singler, sprintf("singler_%s_percell.rds", short)))

    col_name_cell <- sprintf("singler_%s_cell", short)
    sobj@meta.data[[col_name_cell]] <- pred_cell[[ref_name]]$labels[
      match(colnames(sobj), rownames(pred_cell[[ref_name]]))]
  }
}

# ============================================================================
# 5. Consensus annotation table
# ============================================================================
log_msg("Building consensus annotation table...")

ref_names <- names(pred_cluster)
cluster_ids <- rownames(pred_cluster[[ref_names[1]]])

cluster_anno <- data.frame(cluster = cluster_ids, stringsAsFactors = FALSE)

for (ref_name in ref_names) {
  short <- refs[[ref_name]]$short
  pc <- pred_cluster[[ref_name]]
  cluster_anno[[paste0(short, "_label")]]  <- pc$labels
  cluster_anno[[paste0(short, "_score")]]  <- round(apply(pc$scores, 1, max), 3)
  cluster_anno[[paste0(short, "_pruned")]] <- pc$pruned.labels
}

# Cell counts
cluster_sizes <- table(sobj$seurat_clusters)
cluster_anno$n_cells <- as.integer(cluster_sizes[cluster_anno$cluster])

# Agreement flag (if exactly 2 references)
if (length(ref_names) == 2) {
  s1 <- refs[[ref_names[1]]]$short
  s2 <- refs[[ref_names[2]]]$short
  cluster_anno$refs_agree <- cluster_anno[[paste0(s1, "_label")]] ==
                             cluster_anno[[paste0(s2, "_label")]]
}

cluster_anno <- cluster_anno %>% arrange(as.integer(cluster))
fwrite(cluster_anno, file.path(dirs$reports, "singler_cluster_annotation.csv"))

cat("\n=== SingleR Per-Cluster Annotation ===\n")
print(as.data.frame(cluster_anno), row.names = FALSE)

# ---- Per-cell label distribution within clusters ----
if (run_per_cell && length(pred_cell) > 0) {
  first_short <- refs[[ref_names[1]]]$short
  col_cell <- sprintf("singler_%s_cell", first_short)

  if (col_cell %in% colnames(sobj@meta.data)) {
    cell_label_dist <- sobj@meta.data %>%
      group_by(seurat_clusters, .data[[col_cell]]) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(seurat_clusters) %>%
      mutate(pct = round(n / sum(n) * 100, 1)) %>%
      arrange(seurat_clusters, desc(pct))
    fwrite(cell_label_dist, file.path(dirs$reports, "percell_label_distribution.csv"))
  }
}

# ============================================================================
# 6. Plots
# ============================================================================
log_msg("Plotting SingleR UMAPs...")

p_cluster <- DimPlot(sobj, group.by = "seurat_clusters", label = TRUE, repel = TRUE,
                     pt.size = 0.1, label.size = 3) +
  labs(title = sprintf("Clusters (res=%s)", 
       sobj@meta.data$seurat_clusters[1] %>% {attr(., "resolution") %||% "?"})) +
  NoLegend()

# ── UMAP per reference ──
ref_plots <- list()
for (ref_name in ref_names) {
  short <- refs[[ref_name]]$short
  col_cluster <- sprintf("singler_%s_cluster", short)

  tryCatch({
    p <- DimPlot(sobj, group.by = col_cluster, label = TRUE, repel = TRUE,
                 pt.size = 0.1, label.size = 3) +
      labs(title = sprintf("SingleR: %s", ref_name)) +
      theme(legend.text = element_text(size = 7))
    ref_plots[[ref_name]] <- p
  }, error = function(e) {
    log_msg(sprintf("  UMAP for %s failed: %s", ref_name, e$message), "warn")
  })
}

if (length(ref_plots) >= 1) {
  tryCatch({
    # Build combined layout
    if (length(ref_plots) == 1) {
      combined <- p_cluster | ref_plots[[1]]
    } else if (length(ref_plots) == 2) {
      combined <- (p_cluster | ref_plots[[1]]) / (p_cluster | ref_plots[[2]])
    } else {
      combined <- wrap_plots(c(list(p_cluster), ref_plots), ncol = 2)
    }

    pdf(file.path(dirs$plots, "06_singler_umap.pdf"), width = 20, height = 16)
    print(combined)
    dev.off()
    png(file.path(dirs$plots, "06_singler_umap.png"), width = 2000, height = 1600, res = 150)
    print(combined)
    dev.off()
    log_msg("  ✅ 06_singler_umap")
  }, error = function(e) log_msg(sprintf("  UMAP plot failed: %s", e$message), "warn"))
}

# ── Score heatmaps per reference ──
log_msg("Plotting score heatmaps...")
for (ref_name in ref_names) {
  short <- refs[[ref_name]]$short
  pc <- pred_cluster[[ref_name]]

  tryCatch({
    pdf(file.path(dirs$plots, sprintf("07_singler_scores_%s.pdf", short)), width = 14, height = 8)
    plotScoreHeatmap(pc, show.pruned = TRUE,
                     main = sprintf("SingleR Scores: %s (per-cluster)", ref_name))
    dev.off()
    log_msg(sprintf("  ✅ 07_singler_scores_%s", short))
  }, error = function(e) log_msg(sprintf("  Score heatmap %s failed: %s", short, e$message), "warn"))
}

# ── Delta distribution per reference ──
for (ref_name in ref_names) {
  short <- refs[[ref_name]]$short
  pc <- pred_cluster[[ref_name]]

  tryCatch({
    pdf(file.path(dirs$plots, sprintf("08_singler_delta_%s.pdf", short)), width = 14, height = 6)
    plotDeltaDistribution(pc, main = sprintf("Delta Distribution: %s", ref_name))
    dev.off()
    log_msg(sprintf("  ✅ 08_singler_delta_%s", short))
  }, error = function(e) log_msg(sprintf("  Delta plot %s failed: %s", short, e$message), "warn"))
}

# ============================================================================
# 7. Save annotated object
# ============================================================================
log_msg("Saving SingleR-annotated object...")
qsave(sobj, file.path(out_base, "objects", "seurat_singler_annotated.qs"), preset = "fast")
saveRDS(sobj, file.path(out_base, "objects", "seurat_singler_annotated.rds"))
log_msg(sprintf("  ✅ seurat_singler_annotated.qs + .rds: %d cells", ncol(sobj)))

rm(sce); gc()

cat("\n")
cat("==============================================================\n")
cat("   SingleR Annotation Complete                                \n")
cat("==============================================================\n")
cat(sprintf("  Species:     %s\n", ANNO_PARAMS$species))
cat(sprintf("  References:  %s\n", paste(names(refs), collapse = ", ")))
cat(sprintf("  Per-cell:    %s\n", ifelse(run_per_cell, "YES", "NO")))
cat(sprintf("  Clusters:    %d\n", length(cluster_ids)))
cat(sprintf("  Cells:       %d\n", ncol(sobj)))
cat(sprintf("  Results  →   %s\n", dirs$singler))
cat(sprintf("  Reports  →   %s\n", dirs$reports))
cat(sprintf("  Plots    →   %s\n", dirs$plots))
cat("==============================================================\n\n")
