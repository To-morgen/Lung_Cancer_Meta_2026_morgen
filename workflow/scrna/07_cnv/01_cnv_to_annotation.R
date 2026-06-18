#!/usr/bin/env Rscript
# ============================================================================
# 01_cnv_to_annotation.R — Integrate SCEVAN labels into main annotated object
#
# This script runs in the MAIN project renv.
# It does NOT require SCEVAN or infercnv — only reads their output RDS.
#
# Input:  results/scrna/07_cnv/scevan/seurat_with_scevan.rds
#         configs/annotation/celltype_mapping.csv
# Output: results/scrna/06_annotate/objects/seurat_annotated_final.rds
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(dplyr)
  library(data.table)
})

source(here("scripts", "utils", "utils_io.R"))

# ---- Paths ----
scevan_obj   <- scrna_base("07_cnv", "scevan", "seurat_with_scevan.rds")
annot_obj    <- scrna_base("06_annotate", "objects", "seurat_annotated.rds")
mapping_file <- scrna_base("configs", "annotation", "celltype_mapping.csv")
out_dir      <- scrna_base("06_annotate", "objects")
report_dir   <- scrna_base("07_cnv", "reports")

cat("\\n")
cat("==============================================================\\n")
cat("   07_cnv: CNV → Annotation Integration                       \\n")
cat("==============================================================\\n\\n")

# ---- Load ----
log_msg("Loading annotated object...")
sobj <- readRDS(annot_obj)
log_msg(sprintf("  %d cells, columns: %s", ncol(sobj),
                paste(head(colnames(sobj@meta.data), 10), collapse=", ")))

log_msg("Loading SCEVAN object...")
scevan <- readRDS(scevan_obj)

# ---- Transfer SCEVAN labels ----
shared_cells <- intersect(colnames(sobj), colnames(scevan))
log_msg(sprintf("Shared cells: %d / %d", length(shared_cells), ncol(sobj)))

scevan_md <- scevan@meta.data[shared_cells, c("scevan_label", "scevan_subclone"), drop = FALSE]
sobj@meta.data[shared_cells, "scevan_label"]    <- scevan_md$scevan_label
sobj@meta.data[shared_cells, "scevan_subclone"] <- scevan_md$scevan_subclone

missing_cells <- setdiff(colnames(sobj), shared_cells)
if (length(missing_cells) > 0) {
  log_msg(sprintf("  %d cells not in SCEVAN → labeled 'unknown'", length(missing_cells)))
  sobj@meta.data[missing_cells, "scevan_label"] <- "unknown"
}

# ---- Confirm tumor-candidate clusters ----
mapping <- fread(mapping_file)
tumor_candidate_labels <- c("Tumor_putative", "Epi_Tumor")
tumor_putative <- mapping$cluster[mapping$celltype_L1 %in% tumor_candidate_labels]
log_msg(sprintf("Tumor-candidate clusters to confirm: %s", paste(tumor_putative, collapse = ", ")))

confirmation <- sobj@meta.data %>%
  filter(seurat_clusters %in% tumor_putative) %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_total = n(),
    n_scevan_tumor = sum(scevan_label == "tumor", na.rm = TRUE),
    pct_scevan_tumor = round(n_scevan_tumor / n() * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_scevan_tumor))

cat("\\n=== SCEVAN Confirmation of Tumor Candidates ===\\n")
print(as.data.frame(confirmation), row.names = FALSE)

confirmed <- confirmation$seurat_clusters[confirmation$pct_scevan_tumor > 50]
log_msg(sprintf("Confirmed tumor (>50%%): clusters %s", paste(confirmed, collapse = ", ")))

sobj@meta.data$celltype_L1[sobj$seurat_clusters %in% confirmed &
                             sobj$celltype_L1 %in% tumor_candidate_labels] <- "Tumor"

remaining <- setdiff(tumor_putative, confirmed)
if (length(remaining) > 0) {
  log_msg(sprintf("⚠️  Still tumor-candidate: clusters %s", paste(remaining, collapse = ", ")))
}

# ---- Convenience columns ----
sobj$is_tumor <- sobj$celltype_L1 == "Tumor"

# ---- Summary ----
final_summary <- sobj@meta.data %>%
  group_by(celltype_L1) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(desc(n_cells))

cat("\\n=== Final celltype_L1 Distribution ===\\n")
print(as.data.frame(final_summary), row.names = FALSE)
fwrite(final_summary, file.path(report_dir, "cnv_annotation_integration_summary.csv"))

# ---- Save ----
out_file <- file.path(out_dir, "seurat_annotated_final.rds")
saveRDS(sobj, out_file)
log_msg(sprintf("✅ Saved: %s (%d cells)", out_file, ncol(sobj)))

cat("\\n==============================================================\\n")
cat("   All downstream analyses should use:                         \\n")
cat(sprintf("   %s\\n", out_file))
cat("==============================================================\\n\\n")
