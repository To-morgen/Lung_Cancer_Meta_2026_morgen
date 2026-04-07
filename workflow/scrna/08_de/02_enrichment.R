#!/usr/bin/env Rscript
# ============================================================================
# 02_enrichment.R — Gene Set Enrichment (fgsea) on pseudobulk DE results
#
# Input:  Per-celltype DE result CSVs from 01_pseudobulk_de.R
# Output: fgsea results per analysis + combined summary
#
# Gene sets from msigdbr (Hallmark, KEGG, Reactome, GO_BP)
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(fgsea)
  library(msigdbr)
  library(data.table)
  library(dplyr)
  library(yaml)
  library(ggplot2)
})

source(here("scripts", "utils", "utils_io.R"))

cat("\\n")
cat("==============================================================\\n")
cat("   Phase 08 Step 2: Gene Set Enrichment (fgsea)               \\n")
cat("==============================================================\\n\\n")

# ---- Config ----
de_config <- yaml::read_yaml(here("configs", "annotation", "de_contrasts.yaml"))
enr_config <- de_config$enrichment

out_base <- here("results", "scrna", "08_de")
enr_dir  <- file.path(out_base, "enrichment")
plot_dir <- file.path(out_base, "plots", "enrichment")
dir.create(enr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load gene sets ----
log_msg("Loading gene sets from msigdbr...")

all_genesets <- list()
for (gs in enr_config$gene_sets) {
  log_msg(sprintf("  %s (collection=%s)", gs$name, gs$collection))

  args <- list(species = gs$species, collection = gs$collection)
  if (!is.null(gs$subcollection) && gs$subcollection != "null") args$subcollection <- gs$subcollection

  msig <- do.call(msigdbr, args)

  # Convert to fgsea format: list of gene sets
  pathways <- split(msig$gene_symbol, msig$gs_name)
  all_genesets[[gs$name]] <- pathways
  log_msg(sprintf("    → %d gene sets", length(pathways)))
}

# ---- Find DE result files ----
deseq2_dir <- file.path(out_base, "deseq2")
de_files <- list.files(deseq2_dir, pattern = "\\.csv$", full.names = TRUE)
# Exclude the combined file
de_files <- de_files[!grepl("all_de_results_combined", de_files)]

log_msg(sprintf("Found %d DE result files", length(de_files)))

# ---- Run fgsea per DE result × gene set collection ----
all_enr_results <- list()
enr_counter <- 0

for (de_file in de_files) {
  fname <- tools::file_path_sans_ext(basename(de_file))
  log_msg(sprintf("\\n--- Enrichment: %s ---", fname))

  de_res <- fread(de_file)

  # Build ranking vector
  if (!"sign_log_p" %in% colnames(de_res)) {
    de_res$sign_log_p <- sign(de_res$log2FoldChange) * -log10(pmax(de_res$pvalue, 1e-300))
  }

  # Remove NA and duplicates
  de_res <- de_res %>%
    filter(!is.na(sign_log_p), !is.na(gene)) %>%
    distinct(gene, .keep_all = TRUE)

  if (nrow(de_res) < 100) {
    log_msg(sprintf("  Only %d genes → skip", nrow(de_res)), "warn")
    next
  }

  # Named ranking vector
  ranks <- setNames(de_res$sign_log_p, de_res$gene)
  ranks <- sort(ranks, decreasing = TRUE)

  for (gs_name in names(all_genesets)) {
    pathways <- all_genesets[[gs_name]]

    fgsea_res <- tryCatch({
      fgsea(pathways = pathways,
            stats = ranks,
            minSize = enr_config$min_size,
            maxSize = enr_config$max_size,
            nPermSimple = enr_config$nPermSimple)
    }, error = function(e) {
      log_msg(sprintf("  fgsea failed for %s × %s: %s", fname, gs_name, e$message), "warn")
      return(NULL)
    })

    if (is.null(fgsea_res) || nrow(fgsea_res) == 0) next

    # Flatten leadingEdge to string
    fgsea_res$leadingEdge_str <- sapply(fgsea_res$leadingEdge, function(x) paste(x, collapse = ";"))
    fgsea_res$leadingEdge <- NULL

    fgsea_res$de_analysis  <- fname
    fgsea_res$geneset_collection <- gs_name

    # Parse axis/celltype/contrast from filename
    parts <- strsplit(fname, "__")[[1]]
    if (length(parts) >= 3) {
      fgsea_res$axis     <- parts[1]
      fgsea_res$celltype <- parts[2]
      fgsea_res$contrast <- parts[3]
    }

    enr_counter <- enr_counter + 1
    all_enr_results[[enr_counter]] <- fgsea_res

    # Save individual
    enr_fname <- sprintf("%s__%s.csv", fname, gs_name)
    fwrite(fgsea_res, file.path(enr_dir, enr_fname))

    n_sig <- sum(fgsea_res$padj < 0.05, na.rm = TRUE)
    n_up  <- sum(fgsea_res$padj < 0.05 & fgsea_res$NES > 0, na.rm = TRUE)
    n_dn  <- sum(fgsea_res$padj < 0.05 & fgsea_res$NES < 0, na.rm = TRUE)
    log_msg(sprintf("  %s: %d sig pathways (↑%d ↓%d) / %d tested",
                    gs_name, n_sig, n_up, n_dn, nrow(fgsea_res)))
  }
}

# ---- Combined ----
if (length(all_enr_results) > 0) {
  combined_enr <- rbindlist(all_enr_results, fill = TRUE)
  fwrite(combined_enr, file.path(enr_dir, "all_enrichment_combined.csv"))
  log_msg(sprintf("\\nCombined: %d rows across %d analyses", nrow(combined_enr), enr_counter))

  # Summary: top enriched per contrast × collection
  top_summary <- combined_enr %>%
    filter(padj < 0.05) %>%
    group_by(contrast, geneset_collection) %>%
    slice_min(padj, n = 5) %>%
    select(contrast, geneset_collection, pathway, NES, padj, size) %>%
    arrange(contrast, geneset_collection, padj)

  fwrite(top_summary, file.path(out_base, "reports", "enrichment_top5_per_contrast.csv"))

  # Focus pathways
  if (!is.null(enr_config$custom_focus)) {
    focus <- combined_enr %>%
      filter(pathway %in% enr_config$custom_focus) %>%
      select(de_analysis, axis, celltype, contrast, geneset_collection,
             pathway, NES, pval, padj, size) %>%
      arrange(pathway, contrast)

    fwrite(focus, file.path(out_base, "reports", "enrichment_focus_pathways.csv"))
    log_msg(sprintf("Focus pathways: %d rows", nrow(focus)))

    if (nrow(focus) > 0) {
      cat("\\n=== Focus Pathways (padj < 0.05) ===\\n")
      focus_sig <- focus %>% filter(padj < 0.05)
      if (nrow(focus_sig) > 0) {
        print(as.data.frame(focus_sig %>% select(celltype, contrast, pathway, NES, padj)),
              row.names = FALSE)
      } else {
        cat("  No focus pathways reached padj < 0.05\\n")
      }
    }
  }
}

cat("\\n==============================================================\\n")
cat("   Enrichment analysis complete                                \\n")
cat(sprintf("   %d enrichment runs\\n", enr_counter))
cat(sprintf("   Results: %s\\n", enr_dir))
cat("==============================================================\\n\\n")
