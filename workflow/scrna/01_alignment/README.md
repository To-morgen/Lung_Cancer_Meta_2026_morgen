# Workflow: Upstream (Alignment & Raw Processing)

## Scripts

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `01_cellranger_count.sh` | Run Cell Ranger count for scRNA-seq samples | FASTQ files | Filtered/raw matrices, BAM, QC |
| `02_cellranger_qc_summary.sh` | Generate QC summary table from Cell Ranger outputs | Cell Ranger outs | `qc_summary.csv` + console report |
| `03_cellranger_subset_test.sh` | Validate setup with subsampled FASTQs | FASTQ files | Subset Cell Ranger output |

## Shared Functions

- `functions/cellranger_utils.sh` — Reusable functions for Cell Ranger operations

## Usage Examples

```bash
# Subset test (validate before full run)
bash workflow/upstream/03_cellranger_subset_test.sh \
    --sample A1_1 --n-reads 1000000 --cores 8 --mem 32

# Full run (all samples)
nohup bash workflow/upstream/01_cellranger_count.sh \
    --dataset INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1 \
    --auto --cores 32 --mem 128 \
    > logs/cellranger_full.log 2>&1 &

# QC summary
bash workflow/upstream/02_cellranger_qc_summary.sh --output-subdir full --auto
Interface Contract
All upstream scripts:

Source .env.sh for paths (DS_RAW, CELLRANGER_PATH, CELLRANGER_REF)
Read FASTQ from {DS_RAW}/fastq/{sample_id}/
Write output to {DS_RAW}/cellranger_out/{subdir}/{sample_id}/
Return exit code 0 on success, non-zero on failure
Support --help flag
Support --auto for sample auto-detection
Support idempotent re-runs (skip completed samples)
