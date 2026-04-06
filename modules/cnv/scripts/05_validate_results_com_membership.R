#!/usr/bin/env Rscript
# DEBUG/INSPECTION ONLY — not called by pipeline

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(data.table)
  library(yaml)
})

# ------------------------------------------------------------------
# Paths derived from project config (no hardcoded relative paths)
# ------------------------------------------------------------------

module_root <- here::here()   # modules/cnv/
proj_root   <- Sys.getenv("LUNGMETA_ROOT")
if (proj_root == "") proj_root <- normalizePath(file.path(module_root, "..", ".."))

cfg <- yaml::read_yaml(file.path(module_root, "configs", "cnv_params.yaml"))

seurat_path        <- file.path(proj_root, cfg$input$seurat_object)
rdata_path         <- file.path(module_root, cfg$output$module_results,
                                "scevan", "output",
                                "LLC_tumor_CNAmtxSubclones.RData")
reference_clusters <- as.character(cfg$scevan$reference_clusters)
cluster_col        <- if (!is.null(cfg$input$cluster_column))
                        cfg$input$cluster_column else "seurat_clusters"

# ------------------------------------------------------------------
# Load Seurat
# ------------------------------------------------------------------

sobj <- readRDS(seurat_path)

if (!cluster_col %in% colnames(sobj@meta.data)) {
  stop("Cluster column not found in Seurat metadata: ", cluster_col)
}

all_barcodes <- colnames(sobj)
cluster_vec <- as.character(sobj[[cluster_col]][, 1])

# ------------------------------------------------------------------
# Load RData and extract tumor-candidate barcodes
# ------------------------------------------------------------------

e <- new.env(parent = emptyenv())
nm <- load(rdata_path, envir = e)

if (!"results.com" %in% nm) {
  stop("results.com not found in: ", rdata_path, "\nObjects: ", paste(nm, collapse = ", "))
}

obj <- e[["results.com"]]

if (is.null(colnames(obj))) {
  stop("results.com has no colnames; cannot map barcodes.")
}

tumor_candidate_barcodes <- colnames(obj)
tumor_candidate_barcodes <- unique(tumor_candidate_barcodes[tumor_candidate_barcodes %in% all_barcodes])

cat("Total Seurat cells: ", length(all_barcodes), "\n", sep = "")
cat("Tumor-candidate cells from results.com: ", length(tumor_candidate_barcodes), "\n", sep = "")

# ------------------------------------------------------------------
# Build validation table
# ------------------------------------------------------------------

dt <- data.table(
  barcode = all_barcodes,
  cluster = cluster_vec
)

dt[, in_results_com := barcode %in% tumor_candidate_barcodes]
dt[, is_reference_cluster := cluster %in% reference_clusters]

# ------------------------------------------------------------------
# Summaries
# ------------------------------------------------------------------

cluster_summary <- dt[, .(
  n_total = .N,
  n_in_results_com = sum(in_results_com),
  pct_in_results_com = round(100 * sum(in_results_com) / .N, 1)
), by = cluster][order(suppressWarnings(as.integer(cluster)), cluster)]

reference_summary <- dt[, .(
  n_total = .N,
  n_in_results_com = sum(in_results_com),
  pct_in_results_com = round(100 * sum(in_results_com) / .N, 1)
), by = is_reference_cluster]

cross_tab <- dt[, .N, by = .(cluster, in_results_com)]

out_dir <- file.path(proj_root, cfg$output$main_results, "reports")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(cluster_summary,   file.path(out_dir, "results_com_membership_by_cluster.csv"))
fwrite(reference_summary, file.path(out_dir, "results_com_membership_reference_check.csv"))
fwrite(cross_tab,         file.path(out_dir, "results_com_membership_crosstab.csv"))

cat("\nCluster summary:\n")
print(cluster_summary)

cat("\nReference contamination check:\n")
print(reference_summary)

cat(sprintf("\nSaved to: %s\n", out_dir))
