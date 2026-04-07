# Lung Cancer Multi-omics Atlas 2026

**Analyst:** morgen  
**Status:** scRNA-seq Phase 3 complete | Entering downstream analysis  
**Updated:** 2026-04-07

---

## Project Overview

Multi-omics integration of lung cancer datasets focusing on the following tracks:

* **scRNA-seq:** LLC tumor model (GPRIN2-FL / GPRIN2-A1 / mock), 6 samples.
* **Bulk RNA-seq:** 22 public LUAD cohorts (vendor-curated, harmonization planned).

---

## Repository Structure

```text
.
├── configs/                 # Configuration files (YAML, CSV)
│   ├── annotation/          # Cell type mapping, DE contrasts
│   ├── datasets/            # Dataset metadata (public configs)
│   ├── params/              # Analysis parameters (QC, annotation, etc.)
│   └── registry/            # Dataset registry schema
├── configs_private/         # Private dataset configs (not on GitHub)
├── data -> HPC symlink      # Raw + processed data (not tracked)
├── docs/                    # SOPs, setup guides
├── metadata/                # Registry CSV, onboarding reports
├── modules/                 # Isolated analysis modules (independent renv)
│   ├── cnv/                 # SCEVAN / inferCNV (separate dependencies)
│   └── cellchat/            # CellChat (planned)
├── results/                 # All analysis outputs (not tracked)
├── scripts/                 # Project-wide utilities
│   ├── setup/               # Environment setup
│   └── utils/               # IO, plotting, registry helpers
├── workflow/                # Analysis pipelines
│   ├── scrna/               # scRNA-seq (main pipeline)
│   ├── bulk/                # Bulk RNA-seq
│   └── intake/              # Dataset onboarding
├── renv.lock                # R environment lockfile
└── .renvignore              # Excludes modules/results from renv scan
```

---

## scRNA-seq Pipeline Status

| Phase | Directory | Status |
| :--- | :--- | :--- |
| **01 Alignment** | `workflow/scrna/01_alignment/` | **DONE** |
| **02 QC** | `workflow/scrna/02_qc/` | **DONE** |
| **03 Normalize** | `workflow/scrna/03_normalize/` | **DONE** |
| **04 Integrate** | `workflow/scrna/04_integrate/` | **DONE** |
| **05 Cluster** | `workflow/scrna/05_cluster/` | **DONE** |
| **06 Annotate** | `workflow/scrna/06_annotate/` | **DONE** (Final) |
| **07 CNV** | `modules/cnv/` + `workflow/scrna/07_cnv/` | **DONE** |
| **08 DE** | `workflow/scrna/08_de/` | **NEXT** |
| **09 CellChat** | `workflow/scrna/09_cellchat/` | PLANNED |
| **10 Trajectory** | `workflow/scrna/10_trajectory/` | PLANNED |
| **11 Advanced** | `workflow/scrna/11_advanced/` | FUTURE |

---

## Key Design Decisions

1.  **Module Isolation:** Large-scale dependencies like SCEVAN/inferCNV are isolated in `modules/cnv/` with an independent `renv` to prevent dependency hell.
2.  **Adaptive QC:** Implemented **MAD-based QC** (n_mad=3) per sample, utilizing log space for count metrics to account for biological heterogeneity.
3.  **Config-Driven Architecture:** All gene lists, thresholds, and parameters are externalized in **YAML** files — zero hardcoded biological assumptions in the core scripts.
4.  **Iterative Annotation:** Workflow follows an `06_annotate` → `07_cnv` → `06_annotate` loop, using CNV evidence to refine cluster identity.
5.  **Data Integrity:** `seurat_annotated_final.rds` serves as the single "Source of Truth" for all downstream visualization and reporting.

---

## Environment & Execution

* **R Runtime:** `Rscript` via Terminal (Ubuntu HPC).
* **Package Management:** `renv` for reproducible environment locking.
* **Execution Pattern:** `Rscript -e '...'` or direct script execution.
* **Verification:** Visual inspection of plots/results via Cursor.

---

## Quick Start

```bash
# 1. Clone the repository
git clone <repo_url>
cd Lung_Cancer_Meta_2026_morgen

# 2. Establish data link to HPC storage
ln -s /path/to/Lung_Cancer_2026 data

# 3. Restore the R environment
Rscript -e 'renv::restore()'

# 4. Verify environment status
Rscript -e 'renv::status()'
```

---
