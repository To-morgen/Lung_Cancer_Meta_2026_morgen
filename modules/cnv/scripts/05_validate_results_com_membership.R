#!/usr/bin/env Rscript
# ============================================================================
# 05_validate_results_com_membership.R
# DEBUG/INSPECTION ONLY — not called by pipeline
#
# Validates: results.com membership from Subclones .RData
#   → Which clusters are tumor? What's the ref contamination?
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(data.table)
  library(yaml)
})

module_root <- here::here()
proj_root   <- Sys.getenv("LUNGMETA_ROOT")
if (proj_root == "") proj_root <- normalizePath(file.path(module_root, "..", ".."))

cfg <- yaml::read_yaml(file.path(module_root, "configs", "cnv_params.yaml"))

seurat_path        <- file.path(proj_root, cfg$input$seurat_object)
rdata_path         <- file.path(module_root, cfg$output$module_results,
                                "scevan", "output",
                                "LLC_tumor_CNAmtxSubclones.RData")
reference_clusters <- as.character(cfg$scevan$reference_clusters)
cluster_col        <- if (!is.null(cfg$input$cluster_column)) {
                        cfg$input$cluster_column
                      } else {
                        "seurat_clusters"
                      }

cat("\n╔══════════════════════════════════════════════════╗\n")
cat("║  05: Validate results.com membership             ║\n")
cat("╚══════════════════════════════════════════════════╝\n\n")

# ── Load ──
if (!file.exists(rdata_path)) stop("Subclones .RData not found: ", rdata_path)
if (!file.exists(seurat_path)) stop("Seurat object not found: ", seurat_path)

cat(sprintf("Loading .RData:  %s\n", basename(rdata_path)))
env <- new.env()
load(rdata_path, envir = env)
cat(sprintf("Objects loaded: %s\n", paste(ls(env), collapse = ", ")))

cat(sprintf("Loading Seurat:  %s\n", basename(seurat_path)))
sobj <- readRDS(seurat_path)
cell_barcodes <- colnames(sobj)
cat(sprintf("Seurat cells:    %d\n", length(cell_barcodes)))

# ── Extract results.com ──
if (!"results.com" %in% ls(env)) {
  stop("results.com not found in .RData!")
}

results_com <- env$results.com
cat(sprintf("\nresults.com: %d rows, columns: %s\n",
            nrow(results_com), paste(colnames(results_com), collapse = ", ")))

# results.com barcodes
if (!is.null(rownames(results_com))) {
  rc_barcodes <- rownames(results_com)
} else {
  rc_barcodes <- results_com[[1]]
}
cat(sprintf("results.com barcodes: %d\n", length(rc_barcodes)))

# ── Match to Seurat metadata ──
matched <- rc_barcodes[rc_barcodes %in% cell_barcodes]
unmatched <- rc_barcodes[!rc_barcodes %in% cell_barcodes]
cat(sprintf("Matched to Seurat: %d / %d (unmatched: %d)\n",
            length(matched), length(rc_barcodes), length(unmatched)))

# ── Cluster distribution ──
clusters <- as.character(sobj@meta.data[matched, cluster_col])
cluster_tab <- as.data.table(table(cluster = clusters))
setnames(cluster_tab, "N", "n_in_results_com")

# Total cells per cluster
all_clusters <- as.character(sobj@meta.data[[cluster_col]])
total_tab <- as.data.table(table(cluster = all_clusters))
setnames(total_tab, "N", "n_total")

cluster_summary <- merge(total_tab, cluster_tab, by = "cluster", all.x = TRUE)
cluster_summary[is.na(n_in_results_com), n_in_results_com := 0]
cluster_summary[, pct_tumor := round(n_in_results_com / n_total * 100, 2)]
cluster_summary[, is_reference := cluster %in% reference_clusters]
setorder(cluster_summary, -pct_tumor)

cat("\n=== Cluster Membership ===\n")
print(cluster_summary)

# ── Reference contamination ──
ref_in_tumor <- cluster_summary[is_reference == TRUE & n_in_results_com > 0]
total_ref <- cluster_summary[is_reference == TRUE, sum(n_total)]
total_ref_tumor <- cluster_summary[is_reference == TRUE, sum(n_in_results_com)]
ref_contam_pct <- round(total_ref_tumor / total_ref * 100, 2)

reference_summary <- data.table(
  metric = c("total_reference_cells", "reference_in_results_com", "ref_contamination_pct"),
  value = c(total_ref, total_ref_tumor, ref_contam_pct)
)

cat("\n=== Reference Contamination ===\n")
print(reference_summary)
if (ref_contam_pct > 5) {
  cat("⚠️  WARNING: Reference contamination > 5%!\n")
} else {
  cat(sprintf("✅ Reference contamination = %.2f%% (acceptable)\n", ref_contam_pct))
}

# ── Subclone cross-tab (if available) ──
cross_tab <- data.table()
subcl_cols <- grep("subclone|class", colnames(results_com), ignore.case = TRUE, value = TRUE)
if (length(subcl_cols) > 0) {
  cat(sprintf("\nSubclone columns found: %s\n", paste(subcl_cols, collapse = ", ")))
  for (sc in subcl_cols) {
    tab <- as.data.table(table(value = results_com[[sc]]))
    tab$column <- sc
    cross_tab <- rbind(cross_tab, tab)
  }
  cat("\n=== Subclone Distribution ===\n")
  print(cross_tab)
} else {
  cat("\nNo subclone columns found in results.com\n")
  # Try class column on results.com
  if ("class" %in% colnames(results_com)) {
    cross_tab <- as.data.table(table(class = results_com$class))
    cat("\n=== Class Distribution ===\n")
    print(cross_tab)
  }
}

# ── Save ──
out_dir <- file.path(proj_root, cfg$output$main_results, "reports")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(cluster_summary,   file.path(out_dir, "results_com_membership_by_cluster.csv"))
fwrite(reference_summary, file.path(out_dir, "results_com_membership_reference_check.csv"))
if (nrow(cross_tab) > 0) {
  fwrite(cross_tab, file.path(out_dir, "results_com_membership_crosstab.csv"))
}

cat(sprintf("\nOutputs → %s\n", out_dir))
cat("\n========== Validation complete ==========\n")
