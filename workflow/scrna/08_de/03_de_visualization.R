#!/usr/bin/env Rscript
# ============================================================================
# 03_de_visualization.R — Volcano plots + enrichment dotplots
#
# Generates:
#   01. Volcano plots per axis (top contrasts)
#   02. Enrichment dotplot (focus pathways × contrasts)
#   03. DE gene count barplot across celltypes
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

source(here("scripts", "utils", "utils_io.R"))
source(here("workflow", "scrna", "functions", "io_scrna.R"))

out_base <- scrna_base("08_de")
plot_dir <- file.path(out_base, "plots")
enrichment_plot_dir <- file.path(plot_dir, "enrichment")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrichment_plot_dir, recursive = TRUE, showWarnings = FALSE)

cat("\\n")
cat("==============================================================\\n")
cat("   Phase 08 Step 3: DE Visualization                          \\n")
cat("==============================================================\\n\\n")

celltype_display_order <- function(celltype) {
  ifelse(grepl("^pooled_", celltype), 1L, 0L)
}

# ---- Load data ----
combined_file <- file.path(out_base, "deseq2", "all_de_results_combined.csv")
if (!file.exists(combined_file)) stop("Run 01_pseudobulk_de.R first")
de_all <- fread(combined_file)
log_msg(sprintf("Loaded %d DE rows", nrow(de_all)))

# ============================================================================
# Plot 01: Volcano plots for key analyses
# ============================================================================
log_msg("Generating volcano plots...")

make_volcano <- function(df, title, padj_thresh = 0.05, lfc_thresh = 0.5, top_n = 15) {
  # Read from config, with sensible defaults
  if (is.null(padj_thresh)) padj_thresh <- de_config$thresholds$padj %||% 0.05
  if (is.null(lfc_thresh))  lfc_thresh  <- de_config$thresholds$log2fc %||% 0.5
  if (is.null(top_n))       top_n       <- viz_config$volcano$top_n_label %||% 15

  df <- df %>%
    mutate(
      neg_log10_padj = -log10(pmax(padj, 1e-300)),
      color = case_when(
        is.na(padj) ~ "NS",
        padj >= padj_thresh ~ "NS",
        abs(log2FoldChange) < lfc_thresh ~ "NS",
        log2FoldChange > 0 ~ "Up",
        TRUE ~ "Down"
      )
    )

  top_genes <- df %>%
    filter(color != "NS") %>%
    slice_min(padj, n = top_n)

  ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj, color = color)) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70")) +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "grey50") +
    geom_text_repel(data = top_genes, aes(label = gene),
                    size = 2.5, max.overlaps = 20, color = "black") +
    labs(title = title, x = "log2 Fold Change", y = "-log10(padj)", color = "") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Key analyses to plot
key_analyses <- de_all %>%
  distinct(axis, celltype, contrast) %>%
  filter(
    (axis == "tumor_intrinsic" & grepl("pooled|Twist1|Postn", celltype)) |
    (axis == "tme_per_celltype" & grepl("Neutrophil|macrophage|T_cell|NK", celltype, ignore.case = TRUE))
  )

# If too few, just take all pooled + a subset
if (nrow(key_analyses) < 3) {
  key_analyses <- de_all %>%
    distinct(axis, celltype, contrast) %>%
    head(12)
}

volcano_plots <- list()
for (i in seq_len(nrow(key_analyses))) {
  ka <- key_analyses[i, ]
  df <- de_all %>% filter(axis == ka$axis, celltype == ka$celltype, contrast == ka$contrast)
  if (nrow(df) > 0) {
    title <- sprintf("%s | %s\\n%s", ka$celltype, ka$contrast, ka$axis)
    volcano_plots[[paste0(i)]] <- make_volcano(df, title)
  }
}

if (length(volcano_plots) > 0) {
  n_pages <- ceiling(length(volcano_plots) / 6)
  pdf(file.path(plot_dir, "01_volcano_plots.pdf"), width = 18, height = 12)
  for (page in seq_len(n_pages)) {
    idx <- ((page - 1) * 6 + 1):min(page * 6, length(volcano_plots))
    p <- wrap_plots(volcano_plots[idx], ncol = 3)
    print(p)
  }
  dev.off()
  log_msg(sprintf("  Saved: 01_volcano_plots.pdf (%d panels)", length(volcano_plots)))
}

# ============================================================================
# Plot 02: DE gene count barplot
# ============================================================================
log_msg("Generating DEG count barplot...")

summary_file <- file.path(out_base, "reports", "de_summary.csv")
if (file.exists(summary_file)) {
  de_summary <- fread(summary_file)

  de_long <- de_summary %>%
    tidyr::pivot_longer(cols = c(n_up, n_down), names_to = "direction", values_to = "count") %>%
    mutate(
      direction = ifelse(direction == "n_up", "Up", "Down"),
      count_signed = ifelse(direction == "Down", -count, count)
    )

  for (ax in unique(de_long$axis)) {
    ax_data <- de_long %>%
      filter(axis == ax) %>%
      mutate(celltype = factor(
        celltype,
        levels = unique(celltype[order(celltype_display_order(celltype), celltype)])
      ))

    p <- ax_data %>%
      ggplot(aes(x = celltype, y = count_signed, fill = direction)) +
      geom_col() +
      facet_wrap(~contrast, scales = "free_x") +
      scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8")) +
      coord_flip() +
      labs(title = sprintf("DEG counts: %s", ax), x = "", y = "Number of DEGs", fill = "") +
      theme_minimal()

    ggsave(file.path(plot_dir, sprintf("02_deg_counts_%s.pdf", ax)), p, width = 14, height = 10)
    log_msg(sprintf("  Saved: 02_deg_counts_%s.pdf", ax))
  }
}

# ============================================================================
# Plot 03: Enrichment dotplot (focus pathways)
# ============================================================================
log_msg("Generating enrichment dotplot...")

focus_file <- file.path(out_base, "reports", "enrichment_focus_pathways.csv")
if (file.exists(focus_file)) {
  focus <- fread(focus_file)

  if (nrow(focus) > 0) {
    # Clean pathway names
    focus <- focus %>%
      mutate(
        pathway_short = gsub("^HALLMARK_", "", pathway),
        pathway_short = gsub("_", " ", pathway_short),
        sig = ifelse(padj < 0.05, "*", "")
      )

    for (ax in unique(focus$axis)) {
      ax_focus <- focus %>%
        filter(axis == ax) %>%
        mutate(celltype = factor(
          celltype,
          levels = unique(celltype[order(celltype_display_order(celltype), celltype)])
        ))

      p <- ax_focus %>%
        ggplot(aes(x = contrast, y = pathway_short, size = -log10(padj), color = NES)) +
        geom_point() +
        scale_color_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C", midpoint = 0) +
        scale_size_continuous(range = c(1, 6)) +
        facet_wrap(~celltype, scales = "free_x") +
        labs(title = sprintf("Focus Pathways: %s", ax),
             x = "", y = "", color = "NES", size = "-log10(padj)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
              axis.text.y = element_text(size = 7))

      ggsave(file.path(enrichment_plot_dir, sprintf("03_enrichment_focus_%s.pdf", ax)), p, width = 16, height = 12)
      log_msg(sprintf("  Saved: enrichment/03_enrichment_focus_%s.pdf", ax))
    }
  }
} else {
  log_msg("  No focus pathways file → skip", "warn")
}

# ============================================================================
# Plot 04: Enrichment heatmap — top pathways across all tumor contrasts
# ============================================================================
log_msg("Generating enrichment heatmap...")

enr_combined_file <- file.path(out_base, "enrichment", "all_enrichment_combined.csv")
if (file.exists(enr_combined_file)) {
  enr_all <- fread(enr_combined_file)

  # Top Hallmark pathways from tumor analyses
  tumor_hallmark <- enr_all %>%
    filter(grepl("tumor", axis), geneset_collection == "Hallmark", padj < 0.1) %>%
    group_by(pathway) %>%
    summarise(min_padj = min(padj), .groups = "drop") %>%
    slice_min(min_padj, n = 30) %>%
    pull(pathway)

  if (length(tumor_hallmark) > 0) {
    hm_data <- enr_all %>%
      filter(pathway %in% tumor_hallmark, geneset_collection == "Hallmark") %>%
      mutate(label = paste(celltype, contrast, sep = "\\n")) %>%
      select(pathway, celltype, label, NES, padj) %>%
      mutate(label = factor(
        label,
        levels = unique(label[order(celltype_display_order(celltype), celltype, label)])
      ))

    p <- hm_data %>%
      mutate(
        pathway_short = gsub("^HALLMARK_", "", pathway),
        pathway_short = gsub("_", " ", pathway_short),
        sig = ifelse(padj < 0.05, "*", "")
      ) %>%
      ggplot(aes(x = label, y = pathway_short, fill = NES)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = sig), size = 3) +
      scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C", midpoint = 0) +
      labs(title = "Hallmark Pathways: Tumor-intrinsic", x = "", y = "") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            axis.text.y = element_text(size = 7))

    ggsave(file.path(enrichment_plot_dir, "04_hallmark_heatmap_tumor.pdf"), p, width = 16, height = 14)
    log_msg("  Saved: enrichment/04_hallmark_heatmap_tumor.pdf")
  }
}

cat("\\n==============================================================\\n")
cat("   DE Visualization complete                                   \\n")
cat(sprintf("   Plots: %s\\n", plot_dir))
cat(sprintf("   Enrichment plots: %s\\n", enrichment_plot_dir))
cat("==============================================================\\n\\n")
