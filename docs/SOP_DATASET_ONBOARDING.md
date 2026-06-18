下载
→ 注册（registry）
→ intake QC
→ harmonize
→ analysis QC
→ integrate into Atlas
→ run analysis
→ release / archive


# SOP: Dataset Onboarding for Lung Cancer Multi-omics Atlas 2026

## 1. Purpose

This SOP defines the standard workflow for onboarding external datasets into the **Lung Cancer Multi-omics Atlas 2026** project.

The goal is to ensure that every dataset is:

- traceable (**provenance-aware**)
- reproducible (**reproducible**)
- quality-controlled (**QC-gated**)
- harmonized for downstream analysis (**analysis-ready**)
- connected across:
  - **Public GitHub** for code and logic
  - **Private HPC Storage** for heavy data objects
  - **Notion** for project management and decision tracking

---

## 2. Scope

This SOP applies to all dataset types, including but not limited to:

- bulk transcriptomics
- single-cell RNA-seq (scRNA-seq)
- spatial transcriptomics (Spatial)
- whole-exome sequencing (WES)
- proteomics
- clinical annotation tables
- vendor-curated external packages

---

## 3. Dataset Classes

All incoming datasets must be assigned to one of the following classes:

### 3.1 public_original
Public raw or near-raw data directly downloaded from public repositories such as GEO, SRA, ArrayExpress, TCGA, Zenodo, etc.

### 3.2 vendor_curated
Externally curated or packaged datasets prepared by third parties, including `.rds`, `.qs`, `.rda`, processed matrices, merged clinical files, immune infiltration summaries, etc.

### 3.3 internally_curated
Datasets standardized, cleaned, or reconstructed by our own pipeline after intake QC and harmonization.

### 3.4 atlas_ready
Datasets that passed QC and harmonization and are approved for integration into Atlas-level analyses.

### 3.5 release_frozen
Datasets or derived objects frozen for manuscript, thesis, figure generation, or a stable release checkpoint.

---

## 4. Lifecycle States

Each dataset must move through the following lifecycle states:

`discovered -> downloaded -> registered -> intake_qc -> harmonized -> analysis_qc -> analysis_ready -> integrated -> released -> archived`

### State descriptions

- **discovered**: dataset has been identified as potentially relevant
- **downloaded**: files have been obtained and stored in HPC
- **registered**: dataset has been entered into registry
- **intake_qc**: file integrity and structural QC are being checked
- **harmonized**: sample IDs, clinical variables, and object formats are standardized
- **analysis_qc**: post-harmonization QC is being performed
- **analysis_ready**: dataset is approved for downstream analysis
- **integrated**: dataset is connected to the Atlas analysis framework
- **released**: frozen for a stable analysis checkpoint
- **archived**: no longer active, but retained for traceability

---

## 5. Storage Rules

### 5.1 Public GitHub (Code / Logic)

GitHub stores only lightweight and reproducible components:

- scripts
- workflow logic
- configs
- metadata registry
- variable dictionaries
- QC specifications
- documentation
- `renv.lock`

GitHub must **not** store:

- raw sequencing data
- large processed matrices
- heavy object files
- regenerable large figures
- bulky intermediate objects

### 5.2 Private HPC Storage (Data / Objects)

Private HPC storage is the canonical location for:

- downloaded source files
- original vendor packages
- curated standardized objects
- atlas-ready inputs
- large derived outputs
- model cache
- logs

### 5.3 Notion (Management / Decisions)

Notion records:

- dataset registry
- analysis runs
- decisions
- issues / bugs
- progress tracking
- inclusion / exclusion rationale

---

## 6. Directory Principles

Within HPC-backed `data/`, dataset placement should follow lifecycle semantics instead of file extension semantics.

**Recommended structure:**

```text
data/
├── raw/
│   ├── public_original/
│   └── vendor_original/
└── processed/
    ├── curated/
    ├── atlas_ready/
    └── releases/
```

**Notes:**
- `raw/public_original/` stores raw public downloads
- `raw/vendor_original/` stores third-party packaged datasets as-received
- `processed/curated/` stores internally standardized objects
- `processed/atlas_ready/` stores approved analysis inputs
- `processed/releases/` stores frozen analysis versions

---

## 7. Required Identifier

Every dataset must have a unique `dataset_id`.

**Naming convention**
`<Class>_<Source>_<Subtype>_<Modality>_<ShortDesc>_<Date>_vX`

**Example:**
`EXAMPLE_Lung_Bulk_PublicCohort_20260101_v1`

This `dataset_id` must be used consistently across:
- HPC folder names
- registry CSV
- QC reports
- Notion records
- analysis run input field
- result directories
- manuscript supplement naming when relevant

---

## 8. Onboarding Workflow

### 8.1 Step 1 — Discover
**Input**
- Reference from literature, GEO, TCGA, public repository, author website, or external curated source.

**Action**
- Create a preliminary entry in Notion Datasets Registry.

**Output**
- Dataset status = `discovered`

**Required records**
- display name
- source
- disease
- modality
- preliminary URL
- owner
- priority (optional)

### 8.2 Step 2 — Download / Acquire
**Input**
- Verified source URL or vendor package

**Action**
- Download files to private HPC storage

**Output**
Files stored under:
- `data/raw/public_original/<dataset_id>/`
or
- `data/raw/vendor_original/<dataset_id>/`

**Required records**
- download date
- source URL
- file list
- README note if needed

**Minimum checks**
- files exist
- file sizes are non-zero
- storage path is correct

### 8.3 Step 3 — Register
**Input**
- Downloaded dataset files

**Action**
- Register dataset in `metadata/registry/dataset_registry.csv` and in Notion Datasets Registry

**Output**
- Dataset status = `registered`

**Required fields**
- `dataset_id`
- `dataset_type`
- `source_name`
- `modality`
- `subtype`
- `hpc_path`
- `owner`
- `claimed sample size`
- `file format`

### 8.4 Step 4 — Intake QC
**Purpose**
- Determine whether the dataset is structurally usable before harmonization.

**Input**
- Original downloaded files

**Action**
- Run intake QC script and generate machine-readable and human-readable QC outputs

**Output**
- manifest file
- md5 table
- structure summary
- cohort presence map
- intake QC report

**Checks**
- file integrity
- readability
- object structure
- cohort consistency
- sample ID availability
- clinical endpoint availability
- duplicated samples
- obvious missingness or malformed tables

**Status update**
- `intake_qc_status` = `pass` / `warning` / `fail`
- `status` = `intake_qc`

**Failure rule**
- If intake QC fails, dataset must not proceed to harmonization until issues are resolved.

### 8.5 Step 5 — Harmonization
**Purpose**
- Convert heterogeneous external objects into internally standardized analysis objects.

**Input**
- Dataset passing intake QC

**Action**
Standardize:
- sample IDs
- cohort names
- clinical variable names
- survival endpoint variables
- matrix orientation
- feature naming
- metadata schema

**Output**
Stored under:
- `data/processed/curated/<modality>/<dataset_id>/`

**Required outputs**
- curated object
- mapping file
- harmonization log
- aligned sample count summary

**Status update**
- `harmonization_status` = `done`
- `status` = `harmonized`

### 8.6 Step 6 — Analysis QC
**Purpose**
- Check whether the harmonized dataset is statistically safe for downstream modeling.

**Input**
- Curated dataset

**Action**
- Run post-harmonization QC

**Output**
- merged sample count report
- missingness summary
- endpoint consistency summary
- feature distribution summary
- optional exploratory plots

**Checks**
- usable sample number
- endpoint consistency
- duplicated IDs after merge
- missingness burden
- zero variance features
- major outliers
- cohort-specific anomalies

**Status update**
- `analysis_qc_status` = `pass` / `warning` / `fail`

### 8.7 Step 7 — Atlas Integration
**Purpose**
- Decide whether and how the dataset enters the Atlas-level analysis system.

**Input**
- Harmonized dataset with completed QC

**Action**
Assign Atlas role:
- discovery
- validation
- exploratory
- excluded

**Output**
Dataset moved or linked into:
- `data/processed/atlas_ready/<subtype>/<modality>/<dataset_id>/`

**Required records**
- atlas role
- inclusion decision
- rationale
- linked decision entry in Notion

**Status update**
- `status` = `integrated`

### 8.8 Step 8 — Analysis Run Registration
**Purpose**
- Ensure every analysis run is traceable.

**Action**
- For every major run, create an entry in Notion Analysis Runs.

**Required fields**
- `run_id`
- `module`
- `script path`
- `input dataset_id`
- `HPC result path`
- `environment`
- `git commit`
- `release status`

### 8.9 Step 9 — Release / Freeze
**Purpose**
- Freeze stable checkpoints for manuscript, thesis, or reproducible milestone outputs.

**Action**
- save stable outputs
- update `renv.lock`
- commit clean Git state
- record release tag or note
- update Notion

**Output**
Frozen objects and results under:
- `data/processed/releases/`
or
- `results/releases/`

**Status update**
- `status` = `released`

### 8.10 Step 10 — Archive
**Purpose**
- Retain historical traceability while removing inactive datasets from active workflows.

**Action**
- Archive deprecated or excluded datasets without deleting provenance.

**Output**
- Archived status and archived path retained

**Status update**
- `status` = `archived`

---

## 9. QC Gate Rules

A dataset can proceed to the next stage only if the required QC gate is passed.

**Gate 1: intake QC**
- Required to move from `registered` to `harmonized`

**Gate 2: analysis QC**
- Required to move from `harmonized` to `analysis_ready` / `integrated`

Datasets with unresolved structural issues must remain blocked.

---

## 10. Notion Mapping

### 10.1 Datasets Registry
- Used for dataset identity, status, path, source, and inclusion tracking

### 10.2 Analysis Runs
- Used for recording actual analysis executions and links to scripts / outputs

### 10.3 Decisions
- Used for inclusion/exclusion logic, endpoint choices, cohort handling, and major design decisions

### 10.4 Issues / Bugs
- Used for symptoms, root causes, fixes, prevention notes, and linked runs

---

## 11. Required Minimum Artifacts Per Dataset

Every onboarded dataset should have, at minimum:

- registry entry
- HPC storage path
- manifest or file inventory
- intake QC output
- harmonized object or explicit failure note
- status in Notion
- owner
- reproducible script path

---

## 12. Recommended File Outputs

Recommended minimal files for each dataset:

- `<dataset_id>_manifest.csv`
- `<dataset_id>_cohort_presence_map.csv`
- `<dataset_id>_sample_alignment_qc.csv`
- `<dataset_id>_harmonization_log.csv`
- `<dataset_id>_QC_report.html`
- `<dataset_id>_analysis_ready.rds`

---

## 13. Decision Logging Rules

The following events must be recorded in Decisions:

- cohort exclusion
- endpoint redefinition
- sample filtering logic
- switching between OS / PFS / DFS
- using vendor-curated data only for validation
- mapping conflicts that require manual adjudication

---

## 14. Issue Logging Rules

The following events must be recorded in Issues / Bugs:

- unreadable file formats
- package compatibility problems
- sample ID mismatches
- missing endpoints
- duplicated samples
- unexpected cohort count differences
- environment-specific execution failures

---

## 15. Reproducibility Rules

- Use the R version provided by the project’s RStudio Web environment
- Do not mix terminal R and RStudio Web R for formal project runs
- Keep heavy data outside Git
- Use `renv` for R environment locking
- Save stable internal objects in formats with strong long-term compatibility where possible
- Treat vendor outputs as external inputs until validated

---

## 16. Owner Responsibilities

The assigned owner is responsible for:

- source verification
- storage placement
- registry completion
- QC execution
- harmonization review
- Notion updates
- decision logging
- issue logging when needed

---

## 17. Revision History

- **v1.0**: Initial onboarding SOP for Lung Cancer Multi-omics Atlas 2026

***
