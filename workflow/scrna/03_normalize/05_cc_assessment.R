#!/usr/bin/env Rscript
# ============================================================================
# 05_cc_assessment.R — Assess cell cycle effect on PCA
#
# Purpose: Generate diagnostic plots to decide if CC regression is needed.
#
# Decision framework:
#   If CC genes dominate PC1-3     → MUST regress
#   If CC separates clusters       → Should regress (or use CC.Difference)
#   If CC is minor effect          → Don't regress (preserve biology)
#   If studying proliferation      → NEVER regress (it's your signal!)
#
# Output: 6 diagnostic plots → plots/cc_assessment/
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

source(here("workflow", "scrna", "functions", "io_scrna.R"))

# ---- Config ----
out_dirs <- scrna_output_dirs("03_normalize")
cc_plot_dir <- file.path(out_dirs$plots, "cc_assessment")
dir.create(cc_plot_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n╔══════════════════════════════════════════╗\n")
cat("║   Step 05: Cell Cycle Assessment         ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

# ---- Load ----
in_file <- file.path(out_dirs$pca, "merged_pca.rds")
if (!file.exists(in_file)) stop("Input not found: ", in_file)
sobj <- readRDS(in_file)
log_msg(sprintf("Loaded: %d cells × %d genes", ncol(sobj), nrow(sobj)))

# ---- Plot 1: PCA colored by Phase ----
p1 <- DimPlot(sobj, reduction = "pca", group.by = "Phase", dims = c(1, 2),
              pt.size = 0.1, shuffle = TRUE) +
  ggtitle("PC1 vs PC2: Cell Cycle Phase") +
  theme_minimal()

p2 <- DimPlot(sobj, reduction = "pca", group.by = "Phase", dims = c(2, 3),
              pt.size = 0.1, shuffle = TRUE) +
  ggtitle("PC2 vs PC3: Cell Cycle Phase") +
  theme_minimal()

p_phase <- p1 + p2
ggsave(file.path(cc_plot_dir, "01_pca_by_phase.png"), p_phase, width = 14, height = 6, dpi = 150)
log_msg("  01_pca_by_phase.png")

# ---- Plot 2: PCA colored by S.Score and G2M.Score ----
p3 <- FeaturePlot(sobj, reduction = "pca", features = "S.Score", dims = c(1, 2),
                  pt.size = 0.1) +
  scale_color_viridis_c() + ggtitle("PC1 vs PC2: S.Score") + theme_minimal()

p4 <- FeaturePlot(sobj, reduction = "pca", features = "G2M.Score", dims = c(1, 2),
                  pt.size = 0.1) +
  scale_color_viridis_c() + ggtitle("PC1 vs PC2: G2M.Score") + theme_minimal()

p_scores <- p3 + p4
ggsave(file.path(cc_plot_dir, "02_pca_cc_scores.png"), p_scores, width = 14, height = 6, dpi = 150)
log_msg("  02_pca_cc_scores.png")

# ---- Plot 3: PCA colored by sample and group ----
p5 <- DimPlot(sobj, reduction = "pca", group.by = "sample_id", dims = c(1, 2),
              pt.size = 0.1, shuffle = TRUE) +
  ggtitle("PC1 vs PC2: Sample") + theme_minimal()

p6 <- DimPlot(sobj, reduction = "pca", group.by = "group", dims = c(1, 2),
              pt.size = 0.1, shuffle = TRUE) +
  ggtitle("PC1 vs PC2: Group") + theme_minimal()

p_batch <- p5 + p6
ggsave(file.path(cc_plot_dir, "03_pca_by_sample_group.png"), p_batch, width = 14, height = 6, dpi = 150)
log_msg("  03_pca_by_sample_group.png")

# ---- Plot 4: Violin of CC scores by sample ----
p7 <- VlnPlot(sobj, features = "S.Score", group.by = "sample_id", pt.size = 0) +
  ggtitle("S.Score by Sample") + NoLegend() + theme_minimal()

p8 <- VlnPlot(sobj, features = "G2M.Score", group.by = "sample_id", pt.size = 0) +
  ggtitle("G2M.Score by Sample") + NoLegend() + theme_minimal()

p_vln <- p7 / p8
ggsave(file.path(cc_plot_dir, "04_cc_scores_violin.png"), p_vln, width = 10, height = 8, dpi = 150)
log_msg("  04_cc_scores_violin.png")

# ---- Plot 5: Phase proportions per sample ----
phase_df <- as.data.frame(table(sobj$sample_id, sobj$Phase))
colnames(phase_df) <- c("sample", "phase", "count")
phase_df$phase <- factor(phase_df$phase, levels = c("G1", "S", "G2M"))

p9 <- ggplot(phase_df, aes(x = sample, y = count, fill = phase)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("G1" = "#2196F3", "S" = "#FF9800", "G2M" = "#F44336")) +
  labs(title = "Cell Cycle Phase Proportions", y = "Proportion", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(cc_plot_dir, "05_phase_proportions.png"), p9, width = 8, height = 5, dpi = 150)
log_msg("  05_phase_proportions.png")

# ---- Plot 6: Check if CC genes load on top PCs ----
loadings <- Loadings(sobj, reduction = "pca")

# Get CC gene names (mouse format)
s_genes   <- stringr::str_to_title(cc.genes.updated.2019$s.genes)
g2m_genes <- stringr::str_to_title(cc.genes.updated.2019$g2m.genes)
cc_genes  <- unique(c(s_genes, g2m_genes))
cc_in_pca <- intersect(cc_genes, rownames(loadings))

# For each PC, calculate mean absolute loading of CC genes vs all genes
cc_loading_summary <- data.frame(
  PC = 1:min(30, ncol(loadings)),
  cc_mean_abs_loading = sapply(1:min(30, ncol(loadings)), function(i) {
    mean(abs(loadings[cc_in_pca, i]))
  }),
  all_mean_abs_loading = sapply(1:min(30, ncol(loadings)), function(i) {
    mean(abs(loadings[, i]))
  })
)
cc_loading_summary$cc_enrichment <- cc_loading_summary$cc_mean_abs_loading / cc_loading_summary$all_mean_abs_loading

p10 <- ggplot(cc_loading_summary, aes(x = PC, y = cc_enrichment)) +
  geom_bar(stat = "identity", fill = "#9C27B0", alpha = 0.7) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  labs(title = "Cell Cycle Gene Enrichment in PCA Loadings",
       subtitle = "Red dashed = 2x enrichment threshold",
       y = "CC loading / mean loading", x = "PC") +
  theme_minimal()

ggsave(file.path(cc_plot_dir, "06_cc_gene_pca_enrichment.png"), p10, width = 8, height = 5, dpi = 150)
log_msg("  06_cc_gene_pca_enrichment.png")

fwrite(cc_loading_summary, file.path(out_dirs$reports, "cc_pca_enrichment.csv"))

# ---- Print decision helper ----
max_enrich <- max(cc_loading_summary$cc_enrichment[1:5])
top_cc_pcs <- cc_loading_summary$PC[cc_loading_summary$cc_enrichment > 2 & cc_loading_summary$PC <= 10]

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║   CELL CYCLE ASSESSMENT SUMMARY                            ║\n")
cat("╠══════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║   Max CC enrichment in PC1-5: %.2f×                        ║\n", max_enrich))
if (length(top_cc_pcs) > 0) {
  cat(sprintf("║   PCs with >2× CC enrichment (1-10): %s              ║\n",
      paste(top_cc_pcs, collapse = ", ")))
  cat("║   ⚠️  CC may be influencing major PCs                      ║\n")
  cat("║   → Consider regressing S.Score + G2M.Score                ║\n")
  cat("║   → Or use CC.Difference to preserve differentiation      ║\n")
} else {
  cat("║   No PCs (1-10) show >2× CC enrichment                    ║\n")
  cat("║   ✅ CC effect is MINOR — no regression recommended        ║\n")
  cat("║   → Proceed to Phase 4 (Harmony) without CC regression    ║\n")
}
cat("╚══════════════════════════════════════════════════════════════╝\n")

cat("\n========== CC Assessment complete ==========\n")
cat(sprintf("All plots → %s\n", cc_plot_dir))

ls_plots <- list.files(cc_plot_dir, pattern = "\\.png$")
cat(sprintf("Generated %d plots:\n", length(ls_plots)))
for (f in ls_plots) cat(sprintf("  %s\n", f))

gc()
