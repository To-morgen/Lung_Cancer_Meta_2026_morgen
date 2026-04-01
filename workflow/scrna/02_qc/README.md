# Phase 2: QC, Ambient RNA Removal, Doublet Detection

## Pipeline Order (updated)

run_qc_pipeline.R
│
├── 00_soupx_ambient.R SoupX ambient RNA removal (on raw CR output)
├── 01_create_seurat.R Create per-sample Seurat objects + QC metrics
├── 02_qc_filter.R MAD-based QC filtering (dead cells, debris)
├── 03_doublet_removal.R scDblFinder on CLEAN cells (after QC)
└── 04_qc_visualization.R 6-panel diagnostic plots


**Key design decision**: Doublet detection runs AFTER QC filtering.
Dead/dying cells (high mito%) confuse doublet simulation.
Ref: scDblFinder documentation recommends minimal QC before detection.

## QC Method: MAD (Median Absolute Deviation)

For each sample independently:
nFeature_RNA: log10 space, ± 3×MAD from median (hard floor: 200)
nCount_RNA: log10 space, ± 3×MAD from median (hard floor: 500)
percent.mt: linear space, +3×MAD from median (hard cap: 20%)


Reference: Lun et al., 2016; Amezquita et al., Orchestrating Single-Cell Analysis

## Output Structure

results/scrna/02_qc/
├── soupx/ {sample}_soupx_corrected.rds
├── raw_seurat/ {sample}_raw.rds (unfiltered, with QC metrics)
├── qc_filtered/ {sample}_qc_filtered.rds (MAD filtered, doublets labeled)
├── doublets/ {sample}_doublet_calls.rds
├── clean/ {sample}_clean.rds (final: QC + doublet removed)
├── plots/ soupx/ + qc/ (PDF + PNG)
└── reports/ soupx_summary.csv, qc_filter_summary.csv,
mad_thresholds.csv, doublet_summary.csv,
raw_seurat_summary.csv


## Configuration Priority

configs_private/datasets/{DATASET_ID}.yaml (dataset-specific)
configs/params/scrna_qc_params.yaml (project defaults)
Hardcoded in functions/qc_utils.R (safety net)
