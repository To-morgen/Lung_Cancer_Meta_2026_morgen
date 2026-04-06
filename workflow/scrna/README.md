# scRNA-seq Workflow

This directory stores the **project-standard scRNA-seq workflow** for the Lung Cancer Multi-omics Atlas 2026 project.

**It is designed for:**
* reproducible staged processing
* stable promotion from raw sample objects to annotated analysis objects
* clean handoff into downstream modules such as CNV, CellChat, and DE analysis

---

## Workflow Stages

### Stage 01 — Alignment (`01_alignment/`)
* Cell Ranger count
* alignment QC summaries
* small subset checks when needed

### Stage 02 — Per-sample QC (`02_qc/`)
* SoupX ambient correction
* Seurat object creation
* MAD-based QC
* doublet detection/removal
* QC visualization

### Stage 03 — Normalization & Feature Extraction (`03_normalize/`)
* sample merge
* cell cycle scoring
* SCTransform
* PCA
* cell cycle effect assessment

### Stage 04 — Integration (`04_integrate/`)
* Harmony integration
* integrated dimensional reduction outputs

### Stage 05 — Clustering (`05_cluster/`)
* neighbor graph
* clustering
* UMAP
* cluster QC and resolution review

### Stage 06 — Annotation (`06_annotate/`)
* marker identification
* SingleR auto-annotation
* manual annotation
* annotated Seurat object generation

### Stage 07 — Figures / Downstream Preparation (`07_figures/`)
* Figure-oriented and reporting-oriented outputs for stable presentation.

### Stage 10 — CNV Handoff
Project-level CNV outputs are promoted to:
```text
results/scrna/10_cnv/
├── infercnv/
├── scevan/
├── consensus/
├── plots/
└── reports/
```

---

## Actual Current State

> **Status:** This repository is no longer at the QC-only stage.

The tree already indicates that integration scripts, clustering scripts, annotation scripts and outputs, and project-level CNV landing directories **already exist**. This workflow should be read as an **active staged pipeline**, not a future scaffold.

---

## Key Biological / Analytical Decisions

| Strategy / Decision | Biological & Analytical Rationale |
| :--- | :--- |
| **Doublet removal before merge** | Doublets are physical artifacts and should be handled per sample whenever possible. |
| **MAD-based QC** | Thresholding should be data-driven rather than fixed by habit alone. |
| **SCTransform over naive log-norm** | Preferred for more stable variance handling in this project. |
| **Score, do not blindly regress, Cell Cycle** | In tumor projects, proliferation is often biology, not nuisance. |
| **Harmony for integration** | Used as a practical PCA-space batch correction strategy. |
| **Multiple clustering resolutions review** | Final annotation should not rest on a single arbitrary clustering parameter. |

---

## Canonical Object Progression

A typical object path follows a strict linear promotion:

`per-sample raw object` -> `QC-clean object` -> `merged object` -> `normalized object` -> `integrated object` -> `clustered object` -> **`annotated object (Canonical Mother Object)`** -> `downstream module input`

---

## Handoff to Modules

### Principle
* `workflow/` owns the standard pipeline.
* `modules/` own specialized downstream analysis.
*(Examples: CNV -> `modules/cnv/`, CellChat -> `modules/cellchat/`)*

### CNV Example
* **Module-native artifacts** remain in: `modules/cnv/results/...`
* **Promoted project outputs** go to: `results/scrna/10_cnv/...`

---

## Function Layers

| Layer | Path | Purpose |
| :---: | :--- | :--- |
| **1** | `scripts/` | Project-wide general utilities |
| **2** | `workflow/scrna/functions/` | scRNA-specific helper functions |
| **3** | `workflow/scrna/<stage>/` | Stage scripts and pipeline entrypoints |

---

## Recommended Practice

- [x] keep **stable logic** in stage scripts
- [x] keep **exploratory code** in `scratch/`
- [x] save **machine-readable outputs** at each major stage
- [x] **log** major runs
- [x] promote **only stable outputs** into project-level result directories

