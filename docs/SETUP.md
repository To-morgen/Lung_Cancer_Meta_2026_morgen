# 🧬 Setup & Governance: HPC Workflow

**Project:** `Lung_Cancer_Meta_2026_morgen`

**Goal:** Maintain a public code/metadata repository while leveraging private, large-scale HPC storage.

---

> [!IMPORTANT]
> ### 🚨 Critical Rule: One R Runtime Only
> 
> 
> **All R code execution in this repo MUST use the R version provided by RStudio Web (Server/Workbench).**
> Do **not** initialize `renv`, install packages, or run analysis from **Terminal R**. Terminal R and RStudio Web often diverge in:
> * **R Binaries** (Versions)
> * **Library Paths**
> * **Bioconductor Compatibility States**
> 
> 
> *This is the single most important environment rule in this project.*

---

## 0️⃣ Storage Philosophy

**Objective:** Reproducible, audit-friendly science without polluting Git with bulky or sensitive artifacts.

### 📊 Git Tracking Policy

| Category | **Track in Git (Commit)** ✅ | **Do NOT Track (Ignore)** ❌ |
| --- | --- | --- |
| **Logic** | `scripts/`, `workflow/`, `configs/` | `logs/`, `tmp/`, core dumps |
| **Metadata** | `metadata/` (Sample sheets, clinical tables) | Raw Sequencing Data (FASTQ/BAM/CRAM) |
| **Environment** | `renv.lock`, `.Rprofile`, `renv/activate.R` | `renv/library/`, cache/staging directories |
| **Objects** | Configuration tables (YAML/JSON) | Large RDS/H5AD objects, VCFs, MTX |
| **IDE** | `README.md`, `docs/`, `*.Rproj` | `.Rproj.user/`, `.Rhistory`, `.RData` |

> [!TIP]
> ### 💡 Rule of Thumb
> 
> 
> * **Git** stores **instructions and receipts** (Code, Metadata, Env, Docs).
> * **HPC** stores **ingredients and heavy outputs** (Raw data, Intermediate files, Results).
> 
> 

---

## 1️⃣ Create the `data` Symlink

All pipelines must interface with the `data/` folder, which acts as a bridge to the private HPC data center.

### 🛠️ Implementation

```bash
# Template
ln -s /path/to/Lung_Cancer_2026 data

# Example (Specific to our HPC)
ln -s /data1/morgen/data_center/Lung_Cancer_2026 data

```

### ✅ Verification

```bash
ls -ld data
# Expected: data -> /data1/morgen/data_center/Lung_Cancer_2026

```

**⚠️ Important Constraints:**

* **No Commits:** Never push raw/processed data to GitHub.
* **Symlink Hygiene:** In `.gitignore`, use `data` (NOT `data/`). Git treats symlinks as files; adding the slash makes Git look inside the folder.

---

## 2️⃣ Recommended Project Structure

```text
Lung_Cancer_Meta_2026_morgen/
├── data/                  # 🔗 Symlink to HPC (not tracked)
├── metadata/              # 📑 Essential small tables
├── configs/               # ⚙️ YAML/JSON/TSV configs
├── workflow/              # 🐍 Snakemake / Pipeline logic
├── scripts/               # 🛠️ Stable, reusable scripts
├── scratch/               # 🧪 Exploratory drafts (clean periodically)
├── results/               # 📊 Analysis outputs (untacked)
├── figures/               # 🎨 Figure exports (untracked)
├── docs/                  # 📝 SOPs, Design docs
├── modules/               # 📦 Isolated analysis modules
│   └── cellchat/
│       ├── cellchat.Rproj
│       ├── scripts/
│       └── ...
├── Lung_Cancer_Meta_2026.Rproj # ⚓ Project Root Anchor
└── README.md

```

---

## 3️⃣ Root Project Environment (Base)

### Why a Root Environment?

It covers the "always-needed" toolbox for:

* Bulk RNA-seq & scRNA-seq basics
* Enrichment analysis & Plotting
* I/O Utilities (`qs2`, `data.table`)

### Initialization (RStudio Web Only)

1. Open `Lung_Cancer_Meta_2026.Rproj`.
2. Run: `source("scripts/00_setup_base_env.R")`.
3. Verify: `renv::status()`.

---

## 4️⃣ Module Environments (Isolation)

### Why Use Modules?

Fragile or dependency-heavy toolchains (e.g., **CellChat**, **SCENIC**, **Arrow**) should be isolated.

* **Avoid** dependency conflicts.
* **Reduce** upgrade risks.
* **Ensure** module-specific reproducibility.

### 🚀 Creating a Module Skeleton

1. **Terminal (Filesystem):**
```bash
cd ~/biohub/projects/Lung_Cancer_Meta_2026_morgen
mkdir -p modules/cellchat/{scripts,notebooks,scratch,results,figures}

```


2. **RStudio Web (Project Setup):**
* Create a new RStudio Project in `modules/cellchat/`.
* Save as `modules/cellchat/cellchat.Rproj`.


3. **RStudio Web (Initialization):**
```r
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::init(bare = TRUE)

```



---

## 5️⃣ Running Scripts: Stable vs. Exploratory

### 🏆 Stable Scripts

* **Location:** `scripts/` or `modules/<name>/scripts/`
* **Naming:** `00_*.R`, `01_*.R` (Pipeline-friendly)
* **Standards:** Clear I/O, documented, and rerunnable.

### 🧪 Exploratory Work

* **Location:** `scratch/` or `modules/<name>/scratch/`
* **Standards:** Messy is OK, but keep it small and clean it periodically.

> [!NOTE]
> **Promotion Rule:** When exploratory code becomes reliable, promote it: `scratch/` ➔ `scripts/`. *This is how chaos slowly becomes civilization.*

---

## 6️⃣ Common Pitfalls & Quick Fixes

### ❌ Problem: `data` still shows in `git status`

* **Cause:** Incorrect `.gitignore` syntax.
* **Fix:** Ensure it says `data`, not `data/`.
* **Check:** `git check-ignore -v data`.

### ❌ Problem: Package not found / Repo mismatch

* **Fix:** Ensure Bioconductor repos are included:
```r
options(repos = BiocManager::repositories())

```



### ❌ Problem: Terminal R and RStudio Web behave differently

* **Fix:** Stop debugging in Terminal R. Always verify `R.version.string` and `.libPaths()` inside **RStudio Web**.

### ❌ Problem: Arrow fails to build

* **Cause:** Outdated system `cmake`.
* **Fix:** Use `conda` or `mamba` to prepare the build toolchain before installing `arrow` within a module.

---

## 📚 Appendices

<details>
<summary><b>Appendix A: Recommended .gitignore</b></summary>

```gitignore
# symlink to HPC data
data

# R / RStudio
.Rproj.user/
.Rhistory
.RData
.Ruserdata

# renv library and staging
renv/library/
renv/staging/

# logs and temporary files
*.log
*.tmp
*.bak
*.swp

# large outputs
results/
figures/

# keep folder structure
!results/.gitkeep
!figures/.gitkeep

```

</details>

<details>
<summary><b>Appendix B: Minimal Onboarding Checklist</b></summary>

1. **Clone** the repository.
2. **Create** the `data` symlink to HPC storage.
3. **Open** `Lung_Cancer_Meta_2026.Rproj` in **RStudio Web**.
4. **Run** the base environment setup (`00_setup_base_env.R`).
5. **Check** `renv::status()`.
6. **Create** module-specific projects only when needed.
7. **Never** run real project analysis from Terminal R.

</details>