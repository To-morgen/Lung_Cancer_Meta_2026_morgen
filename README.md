# Lung Cancer Multi-omics Atlas 2026

Publication-grade multi-omics atlas project for lung cancer, emphasizing:

  * reproducible pipelines
  * environment governance
  * modular analysis
  * analysis-ready outputs
  * public-code / private-data separation

-----

## Overview

This repository is the **public logic and metadata layer** of the project.

**It tracks:**

  * workflow logic
  * analysis scripts
  * configs
  * registry tables
  * documentation
  * environment locks

**It does *not* store:**

  * raw sequencing data
  * bulky intermediate objects
  * large derived matrices
  * sensitive patient-level source data

> **Note:** Heavy data and analysis artifacts live in private HPC-backed storage and are accessed through the project `data` symlink and standardized result directories.

-----

## Project Architecture

This project is intentionally split into **three layers**:

### 1\. `scripts/`

Project-wide reusable utilities and stable helper code.

### 2\. `workflow/`

Project-wide operational workflows for intake, QC, harmonization, promotion, and standard scRNA processing.

### 3\. `modules/`

Isolated scientific analysis modules with their own environments when needed.
*Examples: `modules/cellchat/`, `modules/cnv/`*

**Rule of thumb:**

  * **general reusable logic** → `scripts/`
  * **project-wide standardized pipelines** → `workflow/`
  * **fragile or tool-specific downstream analyses** → `modules/`

-----

## Current Project State

At the time of writing, the repository already contains:

  * dataset governance and onboarding documents
  * root and module-level `renv` environments
  * scRNA workflow through annotation
  * dedicated CNV analysis landing area:
      * `results/scrna/10_cnv/infercnv/`
      * `results/scrna/10_cnv/scevan/`
      * `results/scrna/10_cnv/consensus/`
      * `results/scrna/10_cnv/plots/`
      * `results/scrna/10_cnv/reports/`

This means the repository is no longer just a scaffold; it is an active analysis project with a defined promotion path from module outputs to project-level results.

-----

## Top-Level Structure

```text
.
├── configs/                  # YAML/CSV/TSV configuration files
├── docs/                     # setup, SOPs, governance documents
├── figures/                  # exported figure outputs (untracked/heavy by default)
├── logs/                     # runtime logs
├── metadata/                 # registry tables, reports, lightweight metadata
├── modules/                  # isolated analysis modules
│   ├── cellchat/
│   └── cnv/
├── results/                  # promoted project-level outputs
│   └── scrna/
├── scratch/                  # exploratory drafts
├── scripts/                  # stable reusable project scripts
├── tests/                    # integration/smoke tests
├── workflow/                 # project operational workflows
├── data -> /private/hpc/...  # symlink to HPC-backed storage
├── Lung_Cancer_Meta_2026.Rproj
├── renv.lock
└── README.md
```

-----

## Runtime Contract

### Core Rule

This project uses **one approved R runtime family** for formal runs:

  * package installation
  * `renv` initialization
  * reproducible analysis runs
  * frozen outputs

...must all be tied to the **same R runtime** used by RStudio Web.

### Important Clarification

Terminal execution is allowed. But formal command-line runs **must** use the explicit `Rscript` binary associated with the approved RStudio Web runtime, not an arbitrary shell `Rscript` from another R installation.

**In practice:**

  * ✅ **Allowed:** terminal launch using the approved explicit `Rscript`
  * ❌ **Not Allowed:** mixing unrelated Terminal R and RStudio Web R in the same project lifecycle

-----

## Reproducibility Model

### Root Environment

The root project uses `renv` for the baseline toolbox.

### Module Environments

Dependency-heavy or fragile toolchains should be isolated inside `modules/<name>/`.
*Examples:* CellChat, CNV tools, SCENIC, Arrow-heavy workflows

### Package Installation Preference

**Recommended:**

  * `pak` for installation
  * `renv` for locking and snapshotting

-----

## Data Governance

| Public GitHub Stores | Private HPC Stores | Notion Stores |
| :--- | :--- | :--- |
| • code<br>• configs<br>• metadata<br>• docs<br>• environment lockfiles | • downloaded source files<br>• raw sequencing data<br>• curated large objects<br>• atlas-ready inputs<br>• large derived outputs<br>• model caches<br>• heavy logs | • dataset registry<br>• analysis runs<br>• decisions<br>• issues / bugs<br>• progress tracking |

-----

## Analysis Philosophy

This repository follows a **promotion model**:

  * exploratory code begins in `scratch/`
  * stable logic is promoted into `scripts/` or `workflow/`
  * tool-specific downstream analyses live in `modules/`
  * module outputs are promoted into standardized project-level result directories

*For example, CNV analysis follows:*

  * raw module artifacts remain in `modules/cnv/results/...`
  * promoted project outputs land in `results/scrna/10_cnv/...`

### scRNA Workflow Status

The scRNA branch is organized as a staged workflow:

1.  alignment
2.  QC
3.  normalization
4.  integration
5.  clustering
6.  annotation
7.  figures / downstream

### CNV Integration

The current repository state indicates that annotation outputs already exist and CNV integration has a dedicated project-level landing zone.

### CNV Module Design

The CNV module is intentionally isolated:

  * `modules/cnv/renv.lock`
  * `modules/cnv/scripts/01_scevan.R`
  * `modules/cnv/scripts/02_infercnv.R`

-----

## Output Policy

### Raw artifacts stay in module space

*Examples:*

  * `modules/cnv/results/scevan/output/*.RData`
  * `modules/cnv/results/scevan/output/*.seg`
  * internal heatmaps and tool-native artifacts

### Promoted artifacts go to project space

*Examples:*

  * `results/scrna/10_cnv/scevan/`
  * `results/scrna/10_cnv/reports/`
  * `results/scrna/10_cnv/plots/`

> **Note:** This separation preserves both **tool-native traceability** and **project-level reporting consistency**.

-----

## Quick Start

**1. Clone the repository**

```bash
git clone <your-repo-url>
cd Lung_Cancer_Meta_2026_morgen
```

**2. Create the HPC symlink**

```bash
ln -s /data1/morgen/data_center/Lung_Cancer_2026 data
```

**3. Open the root project once in RStudio Web**

> Use this to initialize and verify the approved runtime.

**4. Initialize the root environment**

```r
source("scripts/setup/00_setup_base_env.R")
renv::status()
```

**5. Run stable scripts from terminal using the approved runtime**
*Example pattern:*

```bash
/path/to/approved/Rscript workflow/scrna/06_annotate/run_annotation_pipeline.R
/path/to/approved/Rscript modules/cnv/scripts/01_scevan.R
```

-----

## Recommended Conventions

### Naming

Use semantic, versioned names when possible: `<Proj>_<Desc>_<Date>_vX`

### Stable vs Exploratory

  * **stable** → `scripts/`, `workflow/`, `modules/*/scripts/`
  * **exploratory** → `scratch/`, `modules/*/scratch/`

### Tracked vs Untracked

  * **Keep in Git:** code, configs, metadata, docs, lockfiles
  * **Keep out of Git:** raw data, heavy objects, large results, caches, temporary files

-----

## Known Governance TODOs

  - [ ] keep root `README` synchronized with actual tree state
  - [ ] keep `workflow/scrna/README.md` synchronized with real stage progress
  - [ ] add module-level smoke tests for high-risk modules such as CNV
  - [ ] record major runs in Notion Analysis Runs
  - [ ] record environment-specific failures in Issues / Bugs

-----

## Key Documents

  * `docs/SETUP.md`
  * `docs/SOP_DATASET_ONBOARDING.md`
  * `workflow/README.md`
  * `workflow/scrna/README.md`

## License

See LICENSE.

-----
