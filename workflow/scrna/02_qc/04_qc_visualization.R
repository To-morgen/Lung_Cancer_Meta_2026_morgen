#!/usr/bin/env Rscript
# ============================================================================
# 04_qc_visualization.R — Comprehensive QC diagnostic plots
#
# Plots:
#   01. Pre-filter violin (nFeature, nCount, percent.mt) per sample
#   02. Feature vs UMI scatter (colored by MT%)
#   03. MAD threshold overlay (showing computed boundaries per sample)
#   04. Doublet score distribution
#   05. Before vs After comparison (raw → QC → QC+doublet)
#   06. Cell count waterfall chart
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(data.table)
})

source(here("workflow", "scrna", "functions", "io_scrna.R"))
source(here("workflow", "scrna", "functions", "qc_utils.R"))
source(here("scripts", "utils", "utils_plotting.R"))

# ---- Config ----
SAMPLES   <- load_sample_list()
QC_PARAMS <- load_qc_params()
out_dirs  <- scrna_output_dirs("02_qc")
plot_dir  <- file.path(out_dirs$plots, "qc")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n╔══════════════════════════════════════════╗\n")
cat("║   Step 04: QC Visualization               ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

# ---- Load all objects ----
raw_list   <- list()
qc_list    <- list()
clean_list <- list()

for (sid in SAMPLES) {
  f1 <- file.path(out_dirs$base, "raw_seurat", sprintf("%s_raw.rds", sid))
  f2 <- file.path(out_dirs$base, "qc_filtered", sprintf("%s_qc_with_doublet_labels.rds", sid))
  f3 <- file.path(out_dirs$base, "clean", sprintf("%s_clean.rds", sid))
  if (file.exists(f1)) raw_list[[sid]]   <- readRDS(f1)
  if (file.exists(f2)) qc_list[[sid]]    <- readRDS(f2)
  if (file.exists(f3)) clean_list[[sid]]  <- readRDS(f3)
}

# ============================================================================
# Plot 1: Pre-filter violin
# ============================================================================
log_msg("Plot 1: Pre-filter violins...")

if (length(raw_list) > 1) {
  merged_raw <- merge(raw_list[[1]], raw_list[-1], add.cell.ids = names(raw_list))
} else {
  merged_raw <- raw_list[[1]]
}

p1a <- VlnPlot(merged_raw, features = "nFeature_RNA", group.by = "sample_id", pt.size = 0) +
  labs(title = "Genes per Cell (raw)") + theme_project() + NoLegend()
p1b <- VlnPlot(merged_raw, features = "nCount_RNA", group.by = "sample_id", pt.size = 0) +
  scale_y_log10() + labs(title = "UMI per Cell (raw, log10)") + theme_project() + NoLegend()
p1c <- VlnPlot(merged_raw, features = "percent.mt", group.by = "sample_id", pt.size = 0) +
  labs(title = "Mito % (raw)") + theme_project() + NoLegend()
p1d <- VlnPlot(merged_raw, features = "percent.ribo", group.by = "sample_id", pt.size = 0) +
  labs(title = "Ribo % (raw)") + theme_project() + NoLegend()

pdf(file.path(plot_dir, "01_violin_raw.pdf"), width = 14, height = 14)
print((p1a | p1b) / (p1c | p1d))
dev.off()
png(file.path(plot_dir, "01_violin_raw.png"), width = 1400, height = 1400, res = 150)
print((p1a | p1b) / (p1c | p1d))
dev.off()
log_msg("  ✅ 01_violin_raw")
rm(merged_raw); gc()

# ============================================================================
# Plot 2: Feature vs UMI scatter with MT% coloring
# ============================================================================
log_msg("Plot 2: Feature vs UMI scatter...")

scatter_plots <- list()
for (sid in names(raw_list)) {
  md <- raw_list[[sid]]@meta.data
  scatter_plots[[sid]] <- ggplot(md, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
    geom_point(size = 0.2, alpha = 0.4) +
    scale_color_viridis_c(limits = c(0, 30), oob = scales::squish) +
    scale_x_log10() +
    labs(title = sid, x = "UMI (log10)", y = "Genes", color = "MT%") +
    theme_project()
}

p2 <- wrap_plots(scatter_plots, ncol = 3) +
  plot_annotation(title = "Feature vs UMI (raw, colored by MT%)",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

pdf(file.path(plot_dir, "02_scatter_raw.pdf"), width = 16, height = 10)
print(p2)
dev.off()
png(file.path(plot_dir, "02_scatter_raw.png"), width = 1600, height = 1000, res = 150)
print(p2)
dev.off()
log_msg("  ✅ 02_scatter_raw")

# ============================================================================
# Plot 3: MAD threshold overlay
# ============================================================================
log_msg("Plot 3: MAD threshold overlay on distributions...")

mad_thresh_file <- file.path(out_dirs$reports, "mad_thresholds.csv")
if (file.exists(mad_thresh_file)) {
  thresh_df <- fread(mad_thresh_file)

  mad_plots <- list()
  for (sid in names(raw_list)) {
    md <- raw_list[[sid]]@meta.data
    sid_thresh <- thresh_df[thresh_df$sample == sid, ]

    # nFeature distribution with MAD boundaries
    nf_thresh <- sid_thresh[sid_thresh$metric == "nFeature_RNA", ]
    mt_thresh <- sid_thresh[sid_thresh$metric == "percent.mt", ]

    p_nf <- ggplot(md, aes(x = nFeature_RNA)) +
      geom_histogram(bins = 100, fill = "steelblue", alpha = 0.6) +
      geom_vline(xintercept = c(nf_thresh$lower_thresh, nf_thresh$upper_thresh),
                 color = "red", linetype = "dashed", linewidth = 0.8) +
      annotate("text", x = nf_thresh$lower_thresh, y = Inf, label = round(nf_thresh$lower_thresh),
               vjust = 2, hjust = -0.1, color = "red", size = 3) +
      annotate("text", x = nf_thresh$upper_thresh, y = Inf, label = round(nf_thresh$upper_thresh),
               vjust = 2, hjust = 1.1, color = "red", size = 3) +
      labs(title = sprintf("%s: nFeature", sid), x = "nFeature_RNA", y = "Count") +
      theme_project()

    p_mt <- ggplot(md, aes(x = percent.mt)) +
      geom_histogram(bins = 100, fill = "coral", alpha = 0.6) +
      geom_vline(xintercept = mt_thresh$upper_thresh,
                 color = "red", linetype = "dashed", linewidth = 0.8) +
      annotate("text", x = mt_thresh$upper_thresh, y = Inf,
               label = sprintf("%.1f%%", mt_thresh$upper_thresh),
               vjust = 2, hjust = -0.1, color = "red", size = 3) +
      labs(title = sprintf("%s: MT%%", sid), x = "percent.mt", y = "Count") +
      theme_project()

    mad_plots[[paste0(sid, "_nf")]] <- p_nf
    mad_plots[[paste0(sid, "_mt")]] <- p_mt
  }

  pdf(file.path(plot_dir, "03_mad_thresholds.pdf"), width = 16, height = 14)
  print(wrap_plots(mad_plots, ncol = 4) +
          plot_annotation(title = "MAD-based QC Thresholds (red dashed lines)",
                          theme = theme(plot.title = element_text(face = "bold", size = 14))))
  dev.off()
  png(file.path(plot_dir, "03_mad_thresholds.png"), width = 1600, height = 1400, res = 150)
  print(wrap_plots(mad_plots, ncol = 4) +
          plot_annotation(title = "MAD-based QC Thresholds (red dashed lines)",
                          theme = theme(plot.title = element_text(face = "bold", size = 14))))
  dev.off()
  log_msg("  ✅ 03_mad_thresholds")
}

# ============================================================================
# Plot 4: Doublet score distribution
# ============================================================================
log_msg("Plot 4: Doublet scores...")

dbl_plots <- list()
for (sid in names(qc_list)) {
  md <- qc_list[[sid]]@meta.data
  if ("scDblFinder_score" %in% colnames(md) && any(!is.na(md$scDblFinder_score))) {
    dbl_plots[[sid]] <- ggplot(md, aes(x = scDblFinder_score, fill = scDblFinder_class)) +
      geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
      scale_fill_manual(values = c("singlet" = "steelblue", "doublet" = "red3")) +
      labs(title = sid, x = "Doublet Score", y = "Count", fill = "") +
      theme_project()
  }
}

if (length(dbl_plots) > 0) {
  pdf(file.path(plot_dir, "04_doublet_scores.pdf"), width = 14, height = 8)
  print(wrap_plots(dbl_plots, ncol = 3) +
          plot_annotation(title = "scDblFinder Doublet Score Distribution",
                          theme = theme(plot.title = element_text(face = "bold", size = 14))))
  dev.off()
  png(file.path(plot_dir, "04_doublet_scores.png"), width = 1400, height = 800, res = 150)
  print(wrap_plots(dbl_plots, ncol = 3) +
          plot_annotation(title = "scDblFinder Doublet Score Distribution",
                          theme = theme(plot.title = element_text(face = "bold", size = 14))))
  dev.off()
}
log_msg("  ✅ 04_doublet_scores")

# ============================================================================
# Plot 5: Three-stage comparison (Raw → QC → Clean)
# ============================================================================
log_msg("Plot 5: Three-stage comparison...")

comp_data <- list()
for (sid in SAMPLES) {
  for (stage_info in list(
    list(lst = raw_list, name = "1_Raw"),
    list(lst = qc_list, name = "2_QC_filtered"),
    list(lst = clean_list, name = "3_Clean")
  )) {
    if (sid %in% names(stage_info$lst)) {
      md <- stage_info$lst[[sid]]@meta.data
      comp_data <- c(comp_data, list(data.frame(
        sample = sid, stage = stage_info$name,
        nFeature = md$nFeature_RNA, nCount = md$nCount_RNA,
        mt = md$percent.mt
      )))
    }
  }
}

if (length(comp_data) > 0) {
  comp_df <- do.call(rbind, comp_data)
  comp_df$stage <- factor(comp_df$stage, levels = c("1_Raw", "2_QC_filtered", "3_Clean"))
  stage_cols <- c("1_Raw" = "#CCCCCC", "2_QC_filtered" = "#FDB462", "3_Clean" = "#4DAF4A")

  pA <- ggplot(comp_df, aes(x = sample, y = nFeature, fill = stage)) +
    geom_violin(scale = "width", alpha = 0.7, draw_quantiles = 0.5) +
    scale_fill_manual(values = stage_cols) +
    labs(title = "Genes per Cell", y = "nFeature_RNA") +
    theme_project() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  pB <- ggplot(comp_df, aes(x = sample, y = nCount, fill = stage)) +
    geom_violin(scale = "width", alpha = 0.7, draw_quantiles = 0.5) +
    scale_fill_manual(values = stage_cols) + scale_y_log10() +
    labs(title = "UMI per Cell (log10)", y = "nCount_RNA") +
    theme_project() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  pC <- ggplot(comp_df, aes(x = sample, y = mt, fill = stage)) +
    geom_violin(scale = "width", alpha = 0.7, draw_quantiles = 0.5) +
    scale_fill_manual(values = stage_cols) +
    labs(title = "Mito %", y = "percent.mt") +
    theme_project() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  pdf(file.path(plot_dir, "05_three_stage_comparison.pdf"), width = 14, height = 14)
  print(pA / pB / pC + plot_layout(guides = "collect"))
  dev.off()
  png(file.path(plot_dir, "05_three_stage_comparison.png"), width = 1400, height = 1400, res = 150)
  print(pA / pB / pC + plot_layout(guides = "collect"))
  dev.off()
}
log_msg("  ✅ 05_three_stage_comparison")

# ============================================================================
# Plot 6: Cell count waterfall
# ============================================================================
log_msg("Plot 6: Cell count waterfall...")

waterfall_data <- data.frame()
for (sid in SAMPLES) {
  n_raw   <- if (sid %in% names(raw_list)) ncol(raw_list[[sid]]) else NA
  n_qc    <- if (sid %in% names(qc_list)) ncol(qc_list[[sid]]) else NA
  n_clean <- if (sid %in% names(clean_list)) ncol(clean_list[[sid]]) else NA
  waterfall_data <- rbind(waterfall_data, data.frame(
    sample = rep(sid, 3),
    stage = c("1_Raw", "2_QC_filtered", "3_Clean"),
    cells = c(n_raw, n_qc, n_clean)
  ))
}

if (nrow(waterfall_data) > 0) {
  waterfall_data$stage <- factor(waterfall_data$stage,
                                  levels = c("1_Raw", "2_QC_filtered", "3_Clean"))

  p6 <- ggplot(waterfall_data, aes(x = sample, y = cells, fill = stage)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_text(aes(label = scales::comma(cells)),
              position = position_dodge(width = 0.9), vjust = -0.3, size = 2.5) +
    scale_fill_manual(values = c("1_Raw" = "#CCCCCC", "2_QC_filtered" = "#FDB462", "3_Clean" = "#4DAF4A")) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = "Cell Count: Raw → QC → Clean", x = "", y = "Cells", fill = "Stage") +
    theme_project()

  pdf(file.path(plot_dir, "06_cell_count_waterfall.pdf"), width = 12, height = 6)
  print(p6)
  dev.off()
  png(file.path(plot_dir, "06_cell_count_waterfall.png"), width = 1200, height = 600, res = 150)
  print(p6)
  dev.off()
}
log_msg("  ✅ 06_cell_count_waterfall")

# ---- Cleanup ----
rm(raw_list, qc_list, clean_list); gc()

cat("\n========== QC Visualization complete ==========\n")
cat(sprintf("Plots → %s\n", plot_dir))
