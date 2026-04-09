# 🧬 Setup & Governance: HPC Workflow

**Project:** `Lung_Cancer_Meta_2026_morgen`
**Goal:** Reproducible scRNA-seq & multi-omics analysis on HPC, with config-driven pipelines and clear code-data separation.

---

## 0️⃣ Storage Philosophy

> **Git** stores **instructions and receipts** (code, configs, metadata, env locks, docs).
> **HPC** stores **ingredients and heavy outputs** (raw data, intermediate objects, results).

### Git Tracking Policy

| Category | **Track** ✅ | **Ignore** ❌ |
|----------|-------------|--------------|
| Logic | `scripts/`, `workflow/`, `configs/` | `logs/`, `tmp/`, core dumps |
| Metadata | `metadata/` (sample sheets, registry) | Raw sequencing data (FASTQ/BAM) |
| Environment | `renv.lock`, `renv/activate.R`, `.Rprofile` | `renv/library/`, `renv/staging/` |
| Config | YAML/CSV configs | Large RDS/H5AD/MTX objects |
| IDE | `*.Rproj`, `README.md`, `docs/` | `.Rproj.user/`, `.Rhistory`, `.RData` |
| Results | — | `results/`, `figures/` (regenerable) |

---

## 1️⃣ Runtime Environment

### R Execution

All R code in this project runs via **terminal `Rscript`**:

```bash
Rscript workflow/scrna/02_qc/00_soupx_ambient.R 2>&1 | tee logs/phase02_00_soupx.log
```

Results and plots are inspected in **Cursor** (or any local editor with image preview).

### R Version

```bash
Rscript --version
# Confirm: R 4.4.x (matches renv.lock)
```

### renv + pak

This project uses **renv** for reproducibility with **pak** as the underlying install engine:

```r
# In .Rprofile or renv/settings.json:
options(renv.config.pak.enabled = TRUE)
```

Benefits of pak:
- Faster parallel downloads
- Better dependency resolution
- Cleaner error messages

### Snakemake

Pipeline orchestration uses Snakemake in a dedicated conda environment:

```bash
conda activate snakemake_env
# Location: /home/morgen/miniconda3/envs/snakemake_env

# Typical usage (local execution)
cd workflow/datasets/novogene_llc/
snakemake --cores 8 -np        # dry-run
snakemake --cores 8            # execute
snakemake phase02_all --cores 8  # single phase
```

---

## 2️⃣ Data Symlink

All pipelines access data through the `data/` symlink:

```bash
# Create (one-time)
ln -s /data1/morgen/data_center/Lung_Cancer_2026 data

# Verify
ls -ld data
# Expected: data -> /data1/morgen/data_center/Lung_Cancer_2026
```

**Constraints:**
- Never commit raw/processed data to Git
- `.gitignore` uses `data` (NOT `data/`) — Git treats symlinks as files

---

## 3️⃣ Project Structure

```text
Lung_Cancer_Meta_2026_morgen/
├── configs/                    ⚙️  YAML/CSV configs (tracked)
│   ├── annotation/             Cell type mappings, DE contrasts
│   ├── datasets/               Dataset identity YAMLs (public)
│   ├── params/                 QC, annotation, pipeline parameters
│   └── registry/               Registry schema
├── configs_private/            🔒 Sensitive dataset configs (tracked — repo is private)
│   └── datasets/               HPC paths, sample details
├── data/                       🔗 Symlink to HPC (NOT tracked)
├── docs/                       📝 SOPs, setup, context
├── figures/                    🎨 Publication figures (NOT tracked)
│   ├── main/
│   └── supplementary/
├── logs/                       📋 Execution logs (NOT tracked)
├── metadata/                   📑 Registry, reports (tracked)
│   ├── registry/
│   └── reports/
├── modules/                    📦 Isolated analysis modules
│   ├── cellchat/               Independent renv
│   └── cnv/                    Independent renv (SCEVAN)
├── results/                    📊 All outputs (NOT tracked)
│   └── scrna/
│       ├── 01_alignment/
│       ├── 02_qc/
│       ├── 03_normalize/
│       ├── 04_integrate/
│       ├── 05_cluster/
│       ├── 06_annotate/
│       ├── 07_cnv/
│       └── 08_de/
├── scratch/                    🧪 Exploratory drafts (NOT tracked)
├── scripts/                    🛠️ Project-wide utilities (tracked)
│   ├── setup/                  Environment setup scripts
│   └── utils/                  Shared R utility functions
├── tests/                      ✅ Integration tests
├── workflow/                   🐍 Pipeline logic (tracked)
│   ├── bulk/                   Bulk RNA-seq workflows
│   ├── datasets/               Per-dataset Snakemake entry points
│   │   ├── novogene_llc/       Snakefile + config.yaml
│   │   └── gse253718_luad/     Snakefile + config.yaml
│   ├── intake/                 Dataset onboarding
│   ├── scrna/                  scRNA-seq phase scripts
│   │   ├── 02_qc/ … 08_de/    Phase scripts (R)
│   │   └── functions/          Shared scRNA R functions
│   └── snakemake/              Shared Snakemake rules
│       ├── rules/              *.smk rule files
│       └── profiles/           Execution profiles (local/slurm)
├── Lung_Cancer_Meta_2026.Rproj
├── renv.lock
└── README.md
```

---

## 4️⃣ Function Layers

```text
Layer 1 — Project-wide utilities          scripts/utils/
  utils_io.R           I/O, logging, config loading
  utils_plotting.R     theme_project(), shared aesthetics
  utils_registry.R     Dataset registry operations

Layer 2 — scRNA modality functions        workflow/scrna/functions/
  io_scrna.R           Matrix loaders, path resolvers, output dir helpers
  qc_utils.R           MAD filtering, QC param loading, config chain
  annotation_utils.R   Annotation config loader, project gene analysis

Layer 3 — Step scripts                    workflow/scrna/{phase}/*.R
  Read config → call Layer 2 functions → write outputs
  Zero hardcoded gene names or thresholds
```

### Config Priority Chain

```text
Dataset private YAML  >  Project params YAML  >  Hardcoded defaults in R functions
(highest priority)       (mid priority)           (safety net)
```

R scripts resolve the active dataset via the `DATASET_ID` environment variable:

```bash
DATASET_ID=novogene_llc Rscript workflow/scrna/02_qc/00_soupx_ambient.R
```

---

## 5️⃣ Module Environments (Isolation)

Fragile or dependency-heavy tools get their own renv:

### Creating a Module

```bash
# 1. Filesystem
mkdir -p modules/cellchat/{scripts,notebooks,scratch,results,figures}

# 2. Create .Rproj
cat > modules/cellchat/cellchat.Rproj << 'EOF'
Version: 1.0

RestoreWorkspace: No
SaveWorkspace: No
AlwaysSaveHistory: Default

EnableCodeIndexing: Yes
UseSpacesForTab: Yes
NumSpacesForTab: 2
Encoding: UTF-8

RnwWeave: Sweave
LaTeX: pdfLaTeX

AutoAppendNewline: Yes
StripTrailingWhitespace: Yes
EOF

# 3. Initialize renv (from within the module directory)
cd modules/cellchat
Rscript -e 'renv::init(bare = TRUE, settings = list(pak.enabled = TRUE))'
```

### Current Modules

| Module | Purpose | renv | Status |
|--------|---------|------|--------|
| `modules/cnv/` | SCEVAN CNV inference | Independent | ✅ Complete |
| `modules/cellchat/` | Cell-cell communication | Independent | ⬜ Placeholder |

---

## 6️⃣ Pipeline Overview

```text
Phase   Directory             Content                          Status
──────────────────────────────────────────────────────────────────────
01      01_alignment          Cell Ranger count (bash)         ✅
02      02_qc                 SoupX → Seurat → MAD → Doublet  ✅
03      03_normalize          Merge → CC → SCT → PCA          ✅
04      04_integrate          Harmony                          ✅
05      05_cluster            UMAP + Cluster QC                ✅
06      06_annotate           Markers + SingleR + Manual       ✅
07      07_cnv                SCEVAN → annotation feedback     ✅
08      08_de                 Pseudobulk DE + Enrichment       ✅
09      09_cellchat           (reserved)                       ⬜
10      10_trajectory         (reserved)                       ⬜
11      11_advanced           (reserved)                       ⬜
```

### Key Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| QC method | MAD (n=3, log10 space) | Data-driven, per-sample adaptive |
| QC → Doublet order | QC first, then scDblFinder | Dead cells confuse doublet simulation |
| Normalization | SCTransform | Better variance stabilization than LogNormalize |
| Cell cycle | Score always; regress CC.Difference only if needed | Preserve proliferation signal in tumor study |
| Batch correction | Harmony | Fast, effective, PCA-space |
| Clustering | Leiden (algorithm 4), multi-resolution | More robust than Louvain |
| DE method | Pseudobulk DESeq2 | Controls pseudo-replication; proper statistical model |
| Annotation | SingleR auto + human-reviewed celltype_mapping.csv | Automation + manual quality gate |
| CNV | SCEVAN in isolated module | Dependency conflicts with main renv |
| Orchestration | Snakemake (per-dataset Snakefile + shared rules) | Reproducible, DAG-aware, dataset-independent |

---

## 7️⃣ Datasets

### Primary: Novogene LLC (mouse)

| Field | Value |
|-------|-------|
| ID | `INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1` |
| Species | Mouse (mm10) |
| Samples | A1_1, A1_2, FL_1, FL_2, mc_1, mc_2 |
| Groups | A1 (GPRIN2 mutant), FL (wildtype), mc (empty vector) |
| Platform | 10x Chromium 3' v3 |
| Status | Phase 2–8 complete |

### Validation: GSE253718 (human)

| Field | Value |
|-------|-------|
| ID | `PUB_GEO_LUAD_scRNA_GSE253718_v1` (planned) |
| Species | Human (GRCh38) |
| Samples | 6 (3 treatment-naive, 3 TKI-resistant) |
| Disease | EGFR-mutant LUAD |
| Format | 10x MTX (matrix.mtx.gz + features.tsv.gz + barcodes.tsv.gz) |
| Status | ⬜ Onboarding planned |

---

## 8️⃣ Snakemake Usage

### Architecture

```text
workflow/datasets/{name}/Snakefile    ← per-dataset entry + targets
workflow/datasets/{name}/config.yaml  ← samples, species, paths
workflow/snakemake/rules/*.smk        ← shared rules (phase02–08)
workflow/scrna/{phase}/*.R            ← analysis logic (unchanged)
```

### Running

```bash
conda activate snakemake_env
cd workflow/datasets/novogene_llc/

# Preview
snakemake --cores 8 -np

# Full pipeline
snakemake --cores 8

# Single phase
snakemake phase02_all --cores 8

# Up to specific output
snakemake results/scrna/05_cluster/objects/seurat_clustered.rds --cores 8

# Visualize DAG
snakemake --dag | dot -Tpdf > pipeline_dag.pdf
```

---

## 9️⃣ Git Workflow

### Branch Convention

```text
main                          Stable, reviewed code
feature/{name}                New functionality
fix/{name}                    Bug fixes
refactor/{name}               Restructuring without behavior change
```

### Typical Commit Flow

```bash
git checkout -b feature/phase09-cellchat
# ... develop ...
git add workflow/scrna/09_cellchat/ configs/...
git status                    # verify NO results/ or data/
git commit -m "feat: Phase 09 CellChat pipeline"
git push -u origin feature/phase09-cellchat
git checkout main
git merge feature/phase09-cellchat
git push origin main
```

### Safety Check Before Push

```bash
# Must return empty (or only .gitkeep)
git ls-files results/ | head
git ls-files logs/ | head

# Confirm gitignore is working
git check-ignore -v results/scrna/02_qc/clean/A1_1_clean.rds
```

---

## 🔟 Common Pitfalls & Fixes

### `data` appears in `git status`

```bash
# .gitignore must say `data` not `data/`
git check-ignore -v data
```

### Package not found

```bash
# Ensure BiocManager repos are included
Rscript -e 'options(repos = BiocManager::repositories()); install.packages("xxx")'
# Or with pak:
Rscript -e 'pak::pkg_install("xxx")'
```

### renv out of sync

```bash
Rscript -e 'renv::status()'
Rscript -e 'renv::restore()'    # install missing
Rscript -e 'renv::snapshot()'   # record new installs
```

### Rplots.pdf appearing in root

```bash
# Add to .gitignore
echo "Rplots.pdf" >> .gitignore
rm -f Rplots.pdf
```

### Module vs root package conflict

```text
Root renv → general scRNA packages (Seurat, harmony, etc.)
Module renv → isolated (SCEVAN, CellChat, etc.)
Never cross-reference library paths between them.
```

---

## Appendix A: Onboarding Checklist (New Collaborator)

```text
□  Clone the repository
□  Create the data symlink: ln -s /data1/morgen/data_center/Lung_Cancer_2026 data
□  Restore root renv: Rscript -e 'renv::restore()'
□  Verify: Rscript -e 'library(Seurat); cat(packageVersion("Seurat"), "\n")'
□  Create snakemake env: conda activate snakemake_env (or install)
□  Dry-run: cd workflow/datasets/novogene_llc && snakemake -np --cores 1
□  Read docs/CONTEXT.md for project background
□  Read docs/SOP_DATASET_ONBOARDING.md for adding new datasets
```

## Appendix B: Creating a New .Rproj

```bash
cat > "project_name.Rproj" << 'EOF'
Version: 1.0

RestoreWorkspace: No
SaveWorkspace: No
AlwaysSaveHistory: Default

EnableCodeIndexing: Yes
UseSpacesForTab: Yes
NumSpacesForTab: 2
Encoding: UTF-8

RnwWeave: Sweave
LaTeX: pdfLaTeX

AutoAppendNewline: Yes
StripTrailingWhitespace: Yes
EOF
```
