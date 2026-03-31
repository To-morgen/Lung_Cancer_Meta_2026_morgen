# Upstream scRNA-seq Pipeline

This directory contains the **from-scratch** single-cell RNA-seq processing pipeline,
designed to take raw FASTQ files through to a registered-ready Seurat object.

## Scope

Handles `internally_generated` datasets — data produced by our own lab and sequenced
by commercial vendors (e.g., Novogene).

## Pipeline Steps

| Step | Script | Tool | Input | Output |
|------|--------|------|-------|--------|
| 0 | `00_create_subset.sh` | seqtk | Full FASTQ | Subset FASTQ (testing) |
| 1 | `01_run_cellranger.sh` | Cell Ranger | FASTQ + reference | per-sample outs/ |
| 2 | `02_ambient_rna_removal.R` | SoupX | raw + filtered matrix | corrected Seurat |
| 3 | `03_doublet_detection.R` | scDblFinder | Seurat object | Seurat + doublet labels |
| 4 | `04_per_sample_qc.R` | Seurat | Seurat + labels | filtered Seurat |
| 5 | `05_merge_and_normalize.R` | Seurat v5 | list of Seurat objects | merged object |
| 6 | `06_to_registered.R` | custom | merged object | registered-ready object |

## Interface Design

Each step communicates through:
- **Manifest CSV files** (metadata interface)
- **Standardized Seurat objects** (data interface)

This means individual tools at each layer can be swapped without breaking downstream steps.

## Configuration

- **Public templates**: `configs/upstream/` (committed to Git)
- **Private instances**: `configs_private/datasets/` (git-ignored, contains sample-level details)

## Execution

All R scripts should be run from **RStudio Web** or via the RStudio-associated `Rscript` binary.
Shell scripts (Cell Ranger) can be run directly in terminal via `nohup`.
