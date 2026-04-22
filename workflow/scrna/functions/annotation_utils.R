# ============================================================================
# annotation_utils.R — Annotation config loader + helpers
# Layer: 2 (modality-level, scRNA-specific)
#
# Species-aware: auto-detects species from dataset config
#   mouse → scrna_annotation_params.yaml
#   human → scrna_annotation_params_human.yaml
# Override: set ANNOTATION_CONFIG env var to force a specific file
# ============================================================================
suppressPackageStartupMessages({
  library(here)
  library(yaml)
})

source(here("scripts", "utils", "utils_io.R"))

# ---- Internal: resolve species from dataset config ----
.resolve_species <- function() {
  # Priority 1: SPECIES env var (explicit override)
  sp <- Sys.getenv("SPECIES", unset = "")
  if (sp != "") {
    log_msg(sprintf("Species from SPECIES env var: %s", sp))
    return(tolower(sp))
  }

  # Priority 2: dataset config species field
  tryCatch({
    # io_scrna.R must be sourced before this
    ds_cfg <- load_dataset_config()
    sp <- ds_cfg$species
    if (!is.null(sp) && sp != "") {
      log_msg(sprintf("Species from dataset config: %s", sp))
      return(tolower(sp))
    }
  }, error = function(e) {
    log_msg(sprintf("Could not read dataset config for species: %s", e$message), "warn")
  })

  # Fallback: mouse

  log_msg("Species not specified — defaulting to mouse", "warn")
  "mouse"
}

#' Load annotation parameters from YAML (species-aware)
#'
#' Resolution order:
#'   1. ANNOTATION_CONFIG env var (full path to YAML)
#'   2. Auto-detect from dataset config species → select matching params file
#'   3. Fallback to scrna_annotation_params.yaml (mouse)
#'
#' @return list with all annotation config
load_annotation_params <- function() {
  # ── Priority 1: explicit env var override ──
  explicit_path <- Sys.getenv("ANNOTATION_CONFIG", unset = "")
  if (explicit_path != "" && file.exists(explicit_path)) {
    params <- yaml::read_yaml(explicit_path)
    log_msg(sprintf("Annotation params loaded [override]: %s", basename(explicit_path)))
    log_msg(sprintf("  %d lineage groups, %d project genes",
                    length(params$lineage_markers), length(params$project_genes)))
    return(params)
  }

  # ── Priority 2: species-based auto-selection ──
  species <- .resolve_species()

  path <- switch(species,
    "human" = here("configs", "params", "scrna_annotation_params_human.yaml"),
    "mouse" = here("configs", "params", "scrna_annotation_params.yaml"),
    # Default for any other species
    here("configs", "params", sprintf("scrna_annotation_params_%s.yaml", species))
  )

  if (!file.exists(path)) {
    # Fallback to default mouse config
    log_msg(sprintf("Species-specific config not found: %s — falling back to mouse", basename(path)), "warn")
    path <- here("configs", "params", "scrna_annotation_params.yaml")
  }

  if (!file.exists(path)) stop("Annotation config not found: ", path)

  params <- yaml::read_yaml(path)

  # Verify species consistency
  config_species <- params$species %||% "unknown"
  if (tolower(config_species) != species) {
    log_msg(sprintf("⚠️  Species mismatch: dataset=%s, annotation_config=%s", species, config_species), "warn")
  }

  log_msg(sprintf("Annotation params loaded [%s]: %s", species, basename(path)))
  log_msg(sprintf("  %d lineage groups, %d project genes",
                  length(params$lineage_markers), length(params$project_genes)))
  params
}

#' Get deduplicated lineage marker vector (ordered, unique)
#' @param params list from load_annotation_params()
#' @return character vector
get_lineage_markers <- function(params) {
  all_genes <- unlist(params$lineage_markers, use.names = FALSE)
  unique(all_genes)
}

#' Get lineage marker reference table (gene + group)
#' @param params list from load_annotation_params()
#' @return data.frame with columns: lineage, gene
get_lineage_reference <- function(params) {
  groups <- params$lineage_markers
  data.frame(
    lineage = rep(names(groups), sapply(groups, length)),
    gene = unlist(groups, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

#' Filter gene list to those present in Seurat object
#' @param genes character vector
#' @param sobj Seurat object
#' @return list(present, missing)
filter_genes_present <- function(genes, sobj) {
  all_genes <- rownames(sobj)
  present <- genes[genes %in% all_genes]
  missing <- setdiff(genes, all_genes)
  if (length(missing) > 0) {
    log_msg(sprintf("  %d/%d genes missing: %s",
                    length(missing), length(genes),
                    paste(head(missing, 10), collapse = ", ")))
  }
  list(present = present, missing = missing)
}

#' Generic per-gene analysis: expression by cluster + group + plots
#' Works for ANY gene list (GPRIN2, EGFR, etc.)
#'
#' @param sobj Seurat object
#' @param genes character vector of gene names
#' @param plot_dir output directory for plots
#' @param report_dir output directory for CSVs
#' @param prefix file name prefix (e.g. "05_project_genes")
analyze_project_genes <- function(sobj, genes, plot_dir, report_dir,
                                  prefix = "05_project_genes") {
  genes_present <- genes[genes %in% rownames(sobj)]
  genes_missing <- setdiff(genes, rownames(sobj))

  if (length(genes_missing) > 0) {
    log_msg(sprintf("  Project genes NOT in data: %s",
                    paste(genes_missing, collapse = ", ")), "warn")
  }
  if (length(genes_present) == 0) {
    log_msg("  No project genes found in data — skipping", "warn")
    return(invisible(NULL))
  }

  md <- sobj@meta.data
  for (gene in genes_present) {
    log_msg(sprintf("  Analyzing project gene: %s", gene))

    expr <- GetAssayData(sobj, assay = "RNA", layer = "counts")[gene, ]
    total_umi <- sum(expr)
    n_pos <- sum(expr > 0)
    n_total <- length(expr)
    log_msg(sprintf("    %.0f total UMI in %.0f/%d cells (%.2f%%)",
                    total_umi, n_pos, n_total, n_pos / n_total * 100))

    # ── Per-group summary ──
    group_summary <- do.call(rbind, lapply(sort(unique(md$group)), function(grp) {
      cells <- rownames(md[md$group == grp, ])
      e <- expr[cells]
      data.frame(gene = gene, group = grp, total_umi = sum(e),
                 n_pos = sum(e > 0), n_total = length(e),
                 pct_pos = round(sum(e > 0) / length(e) * 100, 2),
                 stringsAsFactors = FALSE)
    }))
    data.table::fwrite(group_summary,
      file.path(report_dir, sprintf("%s_expression_by_group.csv", tolower(gene))))

    # ── Per-cluster summary ──
    cl_summary <- md %>%
      dplyr::mutate(.expr = expr[rownames(md)]) %>%
      dplyr::group_by(seurat_clusters) %>%
      dplyr::summarise(
        gene = gene,
        total_umi = sum(.expr),
        n_pos = sum(.expr > 0),
        n_total = dplyr::n(),
        pct_pos = round(sum(.expr > 0) / dplyr::n() * 100, 2),
        mean_expr = round(mean(.expr), 4),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(total_umi))
    data.table::fwrite(cl_summary,
      file.path(report_dir, sprintf("%s_expression_by_cluster.csv", tolower(gene))))

    # ── Plots ──
    if (total_umi > 0 && n_pos > 0) {
      tryCatch({
        p1 <- VlnPlot(sobj, features = gene, group.by = "seurat_clusters", pt.size = 0) +
          labs(title = sprintf("%s by Cluster (%.0f UMI, %.0f cells)", gene, total_umi, n_pos)) +
          NoLegend()
        p2 <- VlnPlot(sobj, features = gene, group.by = "group", pt.size = 0.05) +
          labs(title = sprintf("%s by Group", gene))
        p3 <- FeaturePlot(sobj, features = gene, order = TRUE,
                          cols = c("lightgrey", "darkred"), pt.size = 0.2) +
          labs(title = sprintf("%s on UMAP", gene))

        fname <- sprintf("%s_%s", prefix, tolower(gene))
        pdf(file.path(plot_dir, paste0(fname, ".pdf")), width = 18, height = 12)
        print(p1 / (p2 | p3))
        dev.off()
        png(file.path(plot_dir, paste0(fname, ".png")), width = 1800, height = 1200, res = 150)
        print(p1 / (p2 | p3))
        dev.off()
        log_msg(sprintf("    Done: %s plots", gene))
      }, error = function(e) {
        log_msg(sprintf("    %s plot failed: %s", gene, e$message), "warn")
      })
    } else {
      log_msg(sprintf("    %s has ZERO expression — plots skipped", gene))
      writeLines(
        c(sprintf("%s Expression Report", gene),
          "========================",
          sprintf("%s is in gene list but has 0 UMI across all %d cells.", gene, n_total),
          "",
          "Possible explanations:",
          "  1. Cell line does not express this gene endogenously",
          "  2. Exogenous transgene not in reference genome",
          "  3. Below detection threshold"),
        file.path(report_dir, sprintf("%s_zero_expression_note.txt", tolower(gene)))
      )
    }
  }
  invisible(NULL)
}

cat("[init] workflow/scrna/functions/annotation_utils.R loaded (species-aware)\n")
