#!/usr/bin/env Rscript
# ============================================================================
# 01_pseudobulk_de.R — Pseudobulk Differential Expression (DESeq2)
#
# Input:  seurat_annotated_final.rds + de_contrasts.yaml
# Output: Per-celltype × per-contrast DESeq2 results
#
# Strategy:
#   1. Aggregate raw counts to pseudobulk (per sample × celltype)
#   2. Filter pseudobulk samples (min cells, min UMI, min genes)
#   3. Run DESeq2 for each eligible contrast × celltype
#   4. Export results + summary
#
# ALL parameters from configs/annotation/de_contrasts.yaml
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(DESeq2)
  library(Matrix)
  library(dplyr)
  library(data.table)
  library(yaml)
})

source(here("scripts", "utils", "utils_io.R"))

# ============================================================================
# Config
# ============================================================================
cat("\\n")
cat("==============================================================\\n")
cat("   Phase 08: Pseudobulk Differential Expression (DESeq2)      \\n")
cat("==============================================================\\n\\n")

de_config <- yaml::read_yaml(here("configs", "annotation", "de_contrasts.yaml"))
log_msg("DE config loaded")

out_base <- here("results", "scrna", "08_de")
dirs <- list(
  pseudobulk  = file.path(out_base, "pseudobulk"),
  deseq2      = file.path(out_base, "deseq2"),
  plots       = file.path(out_base, "plots"),
  reports     = file.path(out_base, "reports")
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Helper functions
# ============================================================================

#' Sanitize celltype name for file paths
sanitize_name <- function(x) {
  x <- gsub("[/ ]+", "_", x)
  x <- gsub("[^A-Za-z0-9_.-]", "", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

#' Aggregate to pseudobulk: sum raw counts per sample × celltype
#' @return list(counts = matrix, coldata = data.frame)
aggregate_pseudobulk <- function(sobj, celltype_col, celltype_val,
                                 group_var, sample_var) {
  # Subset to celltype
  cells <- colnames(sobj)[sobj@meta.data[[celltype_col]] %in% celltype_val]
  if (length(cells) == 0) return(NULL)

  sub <- sobj[, cells]
  md <- sub@meta.data

  # Get raw counts from RNA assay
  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")

  # Unique samples
  samples <- sort(unique(md[[sample_var]]))

  # Aggregate
  pb_counts <- matrix(0, nrow = nrow(counts), ncol = length(samples),
                      dimnames = list(rownames(counts), samples))

  for (s in samples) {
    s_cells <- rownames(md[md[[sample_var]] == s, ])
    if (length(s_cells) == 1) {
      pb_counts[, s] <- as.numeric(counts[, s_cells])
    } else if (length(s_cells) > 1) {
      pb_counts[, s] <- Matrix::rowSums(counts[, s_cells])
    }
  }

  # Build coldata
  coldata <- md %>%
    group_by(across(all_of(c(sample_var, group_var)))) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    as.data.frame()
  rownames(coldata) <- coldata[[sample_var]]
  coldata <- coldata[samples, ]  # align order

  list(counts = pb_counts, coldata = coldata)
}

#' Filter pseudobulk samples by min thresholds
filter_pseudobulk <- function(pb, config) {
  min_cells <- config$pseudobulk$min_cells_per_sample
  min_count <- config$pseudobulk$min_total_count
  min_genes <- config$pseudobulk$min_detected_genes

  keep <- rep(TRUE, ncol(pb$counts))
  names(keep) <- colnames(pb$counts)

  # Min cells
  keep <- keep & (pb$coldata$n_cells >= min_cells)

  # Min total count
  total_counts <- colSums(pb$counts)
  keep <- keep & (total_counts >= min_count)

  # Min detected genes
  detected <- colSums(pb$counts > 0)
  keep <- keep & (detected >= min_genes)

  if (sum(keep) < ncol(pb$counts)) {
    removed <- names(keep)[!keep]
    log_msg(sprintf("  Pseudobulk filter: removed %d samples (%s)",
                    length(removed), paste(removed, collapse = ", ")), "warn")
  }

  list(
    counts  = pb$counts[, keep, drop = FALSE],
    coldata = pb$coldata[keep, , drop = FALSE]
  )
}

#' Run DESeq2 for one contrast
#' @return data.frame of results or NULL
run_deseq2_contrast <- function(pb, group_var, contrast_name, num, denom, config) {

  # Check both groups present with ≥ min_samples
  min_samp <- config$pseudobulk$min_samples_per_group
  grp_counts <- table(pb$coldata[[group_var]])
  n_num  <- grp_counts[num]
  n_denom <- grp_counts[denom]

  if (is.na(n_num) || n_num < min_samp) {
    log_msg(sprintf("    %s: skip — %s has %s samples (need %d)",
                    contrast_name, num, ifelse(is.na(n_num), "0", n_num), min_samp), "warn")
    return(NULL)
  }
  if (is.na(n_denom) || n_denom < min_samp) {
    log_msg(sprintf("    %s: skip — %s has %s samples (need %d)",
                    contrast_name, denom, ifelse(is.na(n_denom), "0", n_denom), min_samp), "warn")
    return(NULL)
  }

  # Subset to just these two groups
  keep_samples <- rownames(pb$coldata[pb$coldata[[group_var]] %in% c(num, denom), ])
  cts  <- pb$counts[, keep_samples, drop = FALSE]
  cd   <- pb$coldata[keep_samples, , drop = FALSE]
  cd[[group_var]] <- factor(cd[[group_var]], levels = c(denom, num))

  # Filter very low genes
  keep_genes <- rowSums(cts >= 5) >= 2
  cts <- cts[keep_genes, , drop = FALSE]

  if (nrow(cts) < 100) {
    log_msg(sprintf("    %s: skip — only %d genes pass filter", contrast_name, nrow(cts)), "warn")
    return(NULL)
  }

  # DESeq2
  formula <- as.formula(paste("~", group_var))
  dds <- DESeqDataSetFromMatrix(countData = round(cts),
                                colData = cd,
                                design = formula)

  dds <- tryCatch({
    DESeq(dds, test = config$deseq2$test,
          fitType = config$deseq2$fitType,
          quiet = TRUE)
  }, error = function(e) {
    log_msg(sprintf("    %s: DESeq2 failed — %s", contrast_name, e$message), "ERROR")
    return(NULL)
  })

  if (is.null(dds)) return(NULL)

  # Extract results
  res <- results(dds,
                 contrast = c(group_var, num, denom),
                 alpha = config$deseq2$alpha,
                 independentFiltering = config$deseq2$independent_filtering)

  # Try shrinkage (apeglm needs coefficient name)
  coef_name <- paste0(group_var, "_", num, "_vs_", denom)
  res_shrunk <- tryCatch({
    lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
  }, error = function(e) {
    log_msg(sprintf("    apeglm shrinkage failed, using normal: %s", e$message), "warn")
    tryCatch(
      lfcShrink(dds, coef = coef_name, type = "normal", quiet = TRUE),
      error = function(e2) res  # fallback to unshrunk
    )
  })

  # Build output table
  out <- as.data.frame(res_shrunk) %>%
    tibble::rownames_to_column("gene") %>%
    mutate(
      contrast    = contrast_name,
      numerator   = num,
      denominator = denom,
      # Ranking metric for GSEA: sign(log2FC) × -log10(pvalue)
      sign_log_p  = sign(log2FoldChange) * -log10(pmax(pvalue, 1e-300))
    ) %>%
    arrange(pvalue)

  # Classification
  padj_thresh <- config$thresholds$padj
  lfc_thresh  <- config$thresholds$log2fc
  lfc_strict  <- config$thresholds$log2fc_strict

  out <- out %>%
    mutate(
      significance = case_when(
        is.na(padj) ~ "NS",
        padj >= padj_thresh ~ "NS",
        abs(log2FoldChange) >= lfc_strict ~ "Strong",
        abs(log2FoldChange) >= lfc_thresh ~ "Significant",
        TRUE ~ "NS"
      ),
      direction = case_when(
        significance == "NS" ~ "NS",
        log2FoldChange > 0 ~ "Up",
        TRUE ~ "Down"
      )
    )

  out
}

# ============================================================================
# Load data
# ============================================================================
log_msg("Loading seurat_annotated_final.rds...")
t_load <- Sys.time()
sobj <- readRDS(here("results", "scrna", "06_annotate", "objects", "seurat_annotated_final.rds"))
log_msg(sprintf("  Loaded: %d cells, %d genes (%.1f min)",
                ncol(sobj), nrow(sobj),
                as.numeric(difftime(Sys.time(), t_load, units = "mins"))))

group_var  <- de_config$group_var
sample_var <- de_config$sample_var
contrasts  <- de_config$contrasts

# ============================================================================
# Process each axis
# ============================================================================

all_results    <- list()
all_summaries  <- list()
result_counter <- 0

for (axis_name in names(de_config$axes)) {
  axis <- de_config$axes[[axis_name]]
  log_msg(sprintf("\\n========== Axis: %s ==========", axis_name))
  log_msg(sprintf("  Description: %s", axis$description))

  # Determine which celltypes to analyze
  L1_include <- axis$celltype_L1_include
  L1_exclude <- axis$celltype_L1_exclude

  eligible_cells <- sobj@meta.data %>%
    filter(celltype_L1 %in% L1_include)
  if (!is.null(L1_exclude)) {
    eligible_cells <- eligible_cells %>% filter(!celltype_L1 %in% L1_exclude)
  }

  log_msg(sprintf("  Eligible cells: %d (L1: %s)",
                  nrow(eligible_cells), paste(L1_include, collapse = ", ")))

  # Get L2 types
  L2_types <- sort(unique(eligible_cells$celltype_L2))

  # ── Run per L2 ──
  if (isTRUE(axis$run_per_L2)) {
    for (ct in L2_types) {
      ct_safe <- sanitize_name(ct)
      log_msg(sprintf("\\n  --- %s: %s ---", axis_name, ct))

      # Aggregate
      pb <- aggregate_pseudobulk(sobj, "celltype_L2", ct, group_var, sample_var)
      if (is.null(pb)) { log_msg("    No cells → skip"); next }

      # Filter
      pb <- filter_pseudobulk(pb, de_config)
      if (ncol(pb$counts) < 4) {
        log_msg(sprintf("    Only %d pseudobulk samples after filter → skip", ncol(pb$counts)), "warn")
        next
      }

      log_msg(sprintf("    Pseudobulk: %d samples, %d genes",
                      ncol(pb$counts), nrow(pb$counts)))
      log_msg(sprintf("    Cells per sample: %s",
                      paste(sprintf("%s=%d", pb$coldata[[sample_var]], pb$coldata$n_cells),
                            collapse = ", ")))

      # Save pseudobulk counts
      fwrite(as.data.frame(pb$counts) %>% tibble::rownames_to_column("gene"),
             file.path(dirs$pseudobulk, sprintf("%s__%s_counts.csv", axis_name, ct_safe)))
      fwrite(pb$coldata,
             file.path(dirs$pseudobulk, sprintf("%s__%s_coldata.csv", axis_name, ct_safe)))

      # Run each contrast
      for (con_name in names(contrasts)) {
        con <- contrasts[[con_name]]
        res <- run_deseq2_contrast(pb, group_var, con_name, con$numerator, con$denominator, de_config)

        if (!is.null(res)) {
          res$axis      <- axis_name
          res$celltype  <- ct

          result_counter <- result_counter + 1
          all_results[[result_counter]] <- res

          # Save individual result
          fname <- sprintf("%s__%s__%s.csv", axis_name, ct_safe, con_name)
          fwrite(res, file.path(dirs$deseq2, fname))

          # Summary
          n_up   <- sum(res$direction == "Up", na.rm = TRUE)
          n_down <- sum(res$direction == "Down", na.rm = TRUE)
          n_strong_up   <- sum(res$significance == "Strong" & res$direction == "Up", na.rm = TRUE)
          n_strong_down <- sum(res$significance == "Strong" & res$direction == "Down", na.rm = TRUE)

          log_msg(sprintf("    %s: %d DEGs (↑%d ↓%d | strong ↑%d ↓%d) / %d tested",
                          con_name, n_up + n_down, n_up, n_down,
                          n_strong_up, n_strong_down, nrow(res)))

          all_summaries[[length(all_summaries) + 1]] <- data.frame(
            axis = axis_name, celltype = ct, contrast = con_name,
            n_genes_tested = nrow(res),
            n_sig = n_up + n_down,
            n_up = n_up, n_down = n_down,
            n_strong_up = n_strong_up, n_strong_down = n_strong_down,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  # ── Run pooled ──
  if (isTRUE(axis$run_pooled)) {
    log_msg(sprintf("\\n  --- %s: POOLED ---", axis_name))

    pb <- aggregate_pseudobulk(sobj, "celltype_L1", L1_include, group_var, sample_var)
    if (is.null(pb)) { log_msg("    No cells → skip"); next }

    pb <- filter_pseudobulk(pb, de_config)
    if (ncol(pb$counts) < 4) {
      log_msg(sprintf("    Only %d pseudobulk samples → skip", ncol(pb$counts)), "warn")
      next
    }

    pooled_name <- paste0("pooled_", paste(L1_include, collapse = "_"))
    pooled_safe <- sanitize_name(pooled_name)

    log_msg(sprintf("    Pseudobulk: %d samples, %d genes",
                    ncol(pb$counts), nrow(pb$counts)))

    fwrite(as.data.frame(pb$counts) %>% tibble::rownames_to_column("gene"),
           file.path(dirs$pseudobulk, sprintf("%s__%s_counts.csv", axis_name, pooled_safe)))
    fwrite(pb$coldata,
           file.path(dirs$pseudobulk, sprintf("%s__%s_coldata.csv", axis_name, pooled_safe)))

    for (con_name in names(contrasts)) {
      con <- contrasts[[con_name]]
      res <- run_deseq2_contrast(pb, group_var, con_name, con$numerator, con$denominator, de_config)

      if (!is.null(res)) {
        res$axis     <- axis_name
        res$celltype <- pooled_name

        result_counter <- result_counter + 1
        all_results[[result_counter]] <- res

        fname <- sprintf("%s__%s__%s.csv", axis_name, pooled_safe, con_name)
        fwrite(res, file.path(dirs$deseq2, fname))

        n_up   <- sum(res$direction == "Up", na.rm = TRUE)
        n_down <- sum(res$direction == "Down", na.rm = TRUE)
        n_strong_up   <- sum(res$significance == "Strong" & res$direction == "Up", na.rm = TRUE)
        n_strong_down <- sum(res$significance == "Strong" & res$direction == "Down", na.rm = TRUE)

        log_msg(sprintf("    %s: %d DEGs (↑%d ↓%d | strong ↑%d ↓%d) / %d tested",
                        con_name, n_up + n_down, n_up, n_down,
                        n_strong_up, n_strong_down, nrow(res)))

        all_summaries[[length(all_summaries) + 1]] <- data.frame(
          axis = axis_name, celltype = pooled_name, contrast = con_name,
          n_genes_tested = nrow(res),
          n_sig = n_up + n_down,
          n_up = n_up, n_down = n_down,
          n_strong_up = n_strong_up, n_strong_down = n_strong_down,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # ── Sensitivity analysis ──
  if (!is.null(axis$sensitivity_include_L1)) {
    sens_L1 <- axis$sensitivity_include_L1
    log_msg(sprintf("\\n  --- %s: SENSITIVITY (L1: %s) ---",
                    axis_name, paste(sens_L1, collapse = "+")))

    pb <- aggregate_pseudobulk(sobj, "celltype_L1", sens_L1, group_var, sample_var)
    if (!is.null(pb)) {
      pb <- filter_pseudobulk(pb, de_config)
      if (ncol(pb$counts) >= 4) {
        sens_name <- paste0("sensitivity_", paste(sens_L1, collapse = "_"))
        sens_safe <- sanitize_name(sens_name)

        fwrite(as.data.frame(pb$counts) %>% tibble::rownames_to_column("gene"),
               file.path(dirs$pseudobulk, sprintf("%s__%s_counts.csv", axis_name, sens_safe)))

        for (con_name in names(contrasts)) {
          con <- contrasts[[con_name]]
          res <- run_deseq2_contrast(pb, group_var, con_name, con$numerator, con$denominator, de_config)

          if (!is.null(res)) {
            res$axis     <- axis_name
            res$celltype <- sens_name

            result_counter <- result_counter + 1
            all_results[[result_counter]] <- res

            fname <- sprintf("%s__%s__%s.csv", axis_name, sens_safe, con_name)
            fwrite(res, file.path(dirs$deseq2, fname))

            n_up   <- sum(res$direction == "Up", na.rm = TRUE)
            n_down <- sum(res$direction == "Down", na.rm = TRUE)
            log_msg(sprintf("    %s: %d DEGs (↑%d ↓%d) / %d tested",
                            con_name, n_up + n_down, n_up, n_down, nrow(res)))

            all_summaries[[length(all_summaries) + 1]] <- data.frame(
              axis = axis_name, celltype = sens_name, contrast = con_name,
              n_genes_tested = nrow(res),
              n_sig = n_up + n_down,
              n_up = n_up, n_down = n_down,
              n_strong_up = sum(res$significance == "Strong" & res$direction == "Up", na.rm = TRUE),
              n_strong_down = sum(res$significance == "Strong" & res$direction == "Down", na.rm = TRUE),
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
}

# ============================================================================
# Save combined results + summary
# ============================================================================
log_msg("\\n========== Saving combined outputs ==========")

if (length(all_results) > 0) {
  combined <- rbindlist(all_results, fill = TRUE)
  fwrite(combined, file.path(dirs$deseq2, "all_de_results_combined.csv"))
  log_msg(sprintf("  Combined: %d rows across %d analyses", nrow(combined), length(all_results)))
}

if (length(all_summaries) > 0) {
  summary_df <- do.call(rbind, all_summaries)
  summary_df <- summary_df %>% arrange(axis, celltype, contrast)
  fwrite(summary_df, file.path(dirs$reports, "de_summary.csv"))

  cat("\\n=== DE Summary ===\\n")
  print(as.data.frame(summary_df), row.names = FALSE)
}

cat("\\n==============================================================\\n")
cat("   Pseudobulk DE complete                                      \\n")
cat(sprintf("   %d analyses run across %d axes\\n", result_counter, length(de_config$axes)))
cat(sprintf("   Results: %s\\n", dirs$deseq2))
cat("==============================================================\\n\\n")
