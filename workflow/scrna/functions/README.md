# scRNA-seq Shared Functions (Layer 2)

These functions are shared across all scRNA-seq workflow steps.

## Files

| File | Purpose | Depends on |
|------|---------|------------|
| `io_scrna.R` | Cell Ranger I/O, matrix loading, path helpers | `scripts/utils/utils_io.R` |
| `qc_utils.R` | QC metrics (mito%, ribo%), adaptive thresholds | `scripts/utils/utils_io.R` |
| `cellranger_utils.sh` | Shell-level Cell Ranger utilities | `.env.sh` |

## Layer Hierarchy


Layer 1: scripts/utils/ ← Project-wide (all modalities)
Layer 2: workflow/scrna/functions/ ← scRNA-specific (this directory)
Layer 3: workflow/scrna/02_qc/*.R ← Step-specific scripts


**Rule: Layer N only calls Layer N-1 or lower. Never upward.**
