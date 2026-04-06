#!/usr/bin/env Rscript
# ============================================================================
# 02_infercnv.R — inferCNV: HMM-based CNV inference + publication-grade heatmap
#
# Environment: modules/cnv/ (独立 renv)
# Input:  主项目 annotated Seurat object
# Output: CNV heatmap + HMM predictions → results/scrna/10_cnv/infercnv/
#
# Usage:
#   cd modules/cnv
#   Rscript scripts/02_infercnv.R
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(infercnv)
  library(data.table)
  library(yaml)
})

# ============================================================================
# 1. Paths & Config
# ============================================================================

proj_root <- Sys.getenv("LUNGMETA_ROOT")
if (proj_root == "") {
  proj_root <- normalizePath(file.path(here::here(), "..", ".."))
  cat(sprintf("LUNGMETA_ROOT inferred: %s\n", proj_root))
}

module_root <- here::here()

cat("\n")
cat("╔══════════════════════════════════════════════════╗\n")
cat("║   CNV Module: inferCNV HMM-based Inference       ║\n")
cat(sprintf("║   Time: %s            ║\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("╚══════════════════════════════════════════════════╝\n\n")

cfg <- yaml::read_yaml(file.path(module_root, "configs", "cnv_params.yaml"))
icfg <- cfg$infercnv

cat(sprintf("Analysis mode: %s\n", icfg$analysis_mode))
cat(sprintf("HMM type:      %s\n", icfg$HMM_type))
cat(sprintf("Threads:       %d\n", icfg$num_threads))

# Output dirs
local_out <- file.path(module_root, "results", "infercnv")
main_out  <- file.path(proj_root, cfg$output$main_results, "infercnv")
plot_out  <- file.path(proj_root, cfg$output$main_results, "plots")
report_out <- file.path(proj_root, cfg$output$main_results, "reports")

for (d in c(local_out, main_out, plot_out, report_out)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================================
# 2. Load Seurat
# ============================================================================

input_path <- file.path(proj_root, cfg$input$seurat_object)
cat(sprintf("Loading: %s\n", basename(input_path)))
sobj <- readRDS(input_path)

tryCatch(sobj <- JoinLayers(sobj),
         error = function(e) cat("  JoinLayers skipped (v4 object)\n"))

cat(sprintf("  %d cells, %d genes\n", ncol(sobj), nrow(sobj)))

# ============================================================================
# 3. Prepare Cell Annotations
# ============================================================================

# inferCNV needs a file: barcode \t celltype
# Use cluster numbers as cell type labels (most robust at this stage)
ref_clusters <- as.character(icfg$ref_clusters)

cell_anno <- data.frame(
  barcode  = colnames(sobj),
  celltype = paste0("cluster_", sobj$seurat_clusters),
  stringsAsFactors = FALSE
)
rownames(cell_anno) <- cell_anno$barcode

# Write annotation file
anno_file <- file.path(local_out, "cell_annotations.tsv")
write.table(cell_anno[, "celltype", drop = FALSE],
            file = anno_file, sep = "\t",
            quote = FALSE, col.names = FALSE)

cat(sprintf("Cell annotations: %d cells\n", nrow(cell_anno)))

# Reference group names (cluster-based)
ref_group_names <- paste0("cluster_", ref_clusters)
cat(sprintf("Reference groups: %s\n", paste(ref_group_names, collapse = ", ")))

# ============================================================================
# 4. Gene Order File
# ============================================================================

gene_order_path <- file.path(module_root, icfg$gene_order_file)
if (!file.exists(gene_order_path)) {
  stop("Gene order file not found: ", gene_order_path,
       "\nRun the biomaRt download step first.")
}
cat(sprintf("Gene order file: %s\n", basename(gene_order_path)))

# ============================================================================
# 5. Extract Counts
# ============================================================================

counts <- GetAssayData(sobj, assay = "RNA", layer = "counts")
cat(sprintf("Count matrix: %d genes × %d cells (sparse: %s)\n",
            nrow(counts), ncol(counts), inherits(counts, "sparseMatrix")))

# ============================================================================
# 6. Create inferCNV Object
# ============================================================================

cat("\nCreating inferCNV object...\n")

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file  = anno_file,
  gene_order_file   = gene_order_path,
  ref_group_names   = ref_group_names
)

cat(sprintf("  inferCNV object: %d genes × %d cells\n",
            nrow(infercnv_obj@expr.data), ncol(infercnv_obj@expr.data)))

# ============================================================================
# 7. Run inferCNV
# ============================================================================

cat("\n")
cat("══════════════════════════════════════════════════\n")
cat("   Running inferCNV\n")
cat(sprintf("   Start: %s\n", format(Sys.time(), "%H:%M:%S")))
cat("   ⚠️  This takes 30-120 minutes for 60k cells\n")
cat("══════════════════════════════════════════════════\n\n")

t0 <- Sys.time()

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff              = icfg$cutoff,
  out_dir             = local_out,
  cluster_by_groups   = isTRUE(icfg$cluster_by_groups),
  denoise             = isTRUE(icfg$denoise),
  noise_logFC         = icfg$noise_filter,
  analysis_mode       = icfg$analysis_mode,
  HMM_type            = icfg$HMM_type,
  num_threads         = icfg$num_threads,
  no_prelim_plot      = FALSE,
  no_plot             = FALSE,
  resume_mode         = TRUE    # 如果中断可以续跑
)

elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
cat(sprintf("\n✅ inferCNV completed in %.1f min\n", elapsed))

# ============================================================================
# 8. Save Results
# ============================================================================

cat("\n========== Saving outputs ==========\n")

# Save inferCNV object
saveRDS(infercnv_obj, file.path(local_out, "infercnv_obj_final.rds"))
saveRDS(infercnv_obj, file.path(main_out, "infercnv_obj_final.rds"))

# Copy plots to main project
infercnv_plots <- list.files(local_out, pattern = "\\.(pdf|png)$",
                              full.names = TRUE, recursive = FALSE)
if (length(infercnv_plots) > 0) {
  # Rename with prefix for clarity
  new_names <- file.path(plot_out, paste0("infercnv_", basename(infercnv_plots)))
  file.copy(infercnv_plots, new_names, overwrite = TRUE)
  cat(sprintf("  %d plot(s) → %s\n", length(infercnv_plots), plot_out))
}

# ============================================================================
# 9. Extract CNV Scores Per Cell
# ============================================================================

cat("\nExtracting per-cell CNV scores...\n")

# inferCNV HMM predictions
hmm_files <- list.files(local_out, pattern = "HMM_CNV_predictions",
                         full.names = TRUE, recursive = TRUE)

if (length(hmm_files) > 0) {
  cat(sprintf("  Found %d HMM prediction file(s)\n", length(hmm_files)))
  for (f in hmm_files) {
    fname <- basename(f)
    file.copy(f, file.path(main_out, fname), overwrite = TRUE)
    file.copy(f, file.path(report_out, fname), overwrite = TRUE)
  }
}

# CNV score matrix (if available)
cnv_score_files <- list.files(local_out, pattern = "expr\\.(dat|infercnv)",
                               full.names = TRUE)
if (length(cnv_score_files) > 0) {
  cat(sprintf("  Found CNV expression files: %s\n",
              paste(basename(cnv_score_files), collapse = ", ")))
}

# ============================================================================
# 10. Summary
# ============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════╗\n")
cat("║                inferCNV Complete                     ║\n")
cat(sprintf("║  Elapsed:  %.1f min                                  ║\n", elapsed))
cat(sprintf("║  Cells:    %d                                       ║\n", ncol(sobj)))
cat("╠══════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Module:   modules/cnv/results/infercnv/            ║\n"))
cat(sprintf("║  Main:     results/scrna/10_cnv/infercnv/           ║\n"))
cat(sprintf("║  Plots:    results/scrna/10_cnv/plots/              ║\n"))
cat("╚══════════════════════════════════════════════════════╝\n")

gc()
