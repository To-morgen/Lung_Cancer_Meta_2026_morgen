# **Project Working Memory: Lung Cancer Multi-omics Atlas 2026**

**Version**: 1.0.0  
**Last Updated**: 2026-04-20  
**Analyst**: morgen  
**Status**: Phase 1-8 complete (INT_Novogene), Phase 1-6 complete (GSE253718)

---

## **1. Project Mission**

### **Core Objective**
Build a **config-driven, reproducible, modular scRNA-seq analysis pipeline** for lung cancer multi-omics studies, with emphasis on:
- Zero hardcoded parameters (all in YAML/CSV)
- Phase-based architecture (01-11)
- Multi-dataset support (mouse/human, internal/public)
- Isolated environments for dependency-heavy tools (CNV, CellChat)
- Publication-quality reproducibility (Level 3: external researcher can replicate)

### **Current Focus**
- **Primary dataset**: INT_Novogene_LLC (mouse LLC model, GPRIN2 wildtype/mutant/control, 6 samples)
- **Validation dataset**: GSE253718 (human LUAD, EGFR-mutant, TKI-naive vs resistant, 6 samples)
- **Analysis stage**: Entering Phase 09 (CellChat) after completing Phase 08 (DE + enrichment)

### **Future Plans**
- Meta-analysis across multiple datasets
- Template-ization for rapid project initialization (post-pipeline maturity)
- Comprehensive downstream modules (CellChat, Trajectory, SCENIC, Velocity, Scoring)

---

## **2. Directory Conventions**

### **Top-level Structure**
```
Lung_Cancer_Meta_2026_morgen/
├── configs/                      # Public configs (analysis parameters)
│   ├── annotation/               # Cell type mappings, DE contrasts
│   ├── datasets/                 # Public dataset metadata (GEO accessions)
│   ├── params/                   # QC, normalization, clustering params
│   └── registry/                 # Dataset registry schema
├── configs_private/              # Private configs (HPC paths, unpublished data)
│   └── datasets/                 # Dataset-specific paths and sample manifests
├── data/                         # Symlink to HPC storage (NOT tracked)
├── docs/                         # SOPs, setup guides, pipeline docs
├── metadata/                     # Dataset registry, onboarding reports
│   ├── registry/                 # dataset_registry.csv
│   └── reports/                  # Onboarding QC reports
├── modules/                      # Isolated analysis modules (independent renv)
│   ├── cnv/                      # SCEVAN/inferCNV (separate renv)
│   └── cellchat/                 # CellChat (planned)
├── results/                      # All analysis outputs (NOT tracked)
│   └── scrna/{ds_prefix}/        # Per-dataset results
│       ├── 02_qc/ → 08_de/       # Phase-specific outputs
│       └── .phase{N}_done        # Sentinel files
├── scripts/                      # Project-wide utilities (Layer 1)
│   ├── setup/                    # Environment setup
│   └── utils/                    # IO, plotting, registry helpers
├── workflow/                     # Analysis pipelines
│   ├── datasets/{name}/          # Per-dataset Snakemake entry points
│   │   ├── Snakefile             # Dataset-specific targets
│   │   └── config.yaml           # Samples, species, paths
│   ├── scrna/                    # scRNA-seq main pipeline
│   │   ├── 01_alignment/ → 11_advanced/  # Phase directories
│   │   ├── functions/            # scRNA-specific utilities (Layer 2)
│   │   └── {phase}/run_*_pipeline.R  # Phase orchestrators
│   └── snakemake/rules/          # Shared Snakemake rules
├── renv.lock                     # Main R environment lockfile
└── README.md                     # Project overview
```

### **Key Principles**
- **Git tracks**: Code, configs, metadata, docs, renv.lock
- **Git ignores**: data/, results/, logs/, .Rproj.user/, *.bak
- **HPC stores**: Raw data, processed objects, analysis outputs
- **Symlink**: `data/` → `/data1/morgen/data_center/Lung_Cancer_2026/`

---

## **3. Main Workflow Entry Points**

### **A. Snakemake (Production, Recommended)**
```bash
# Navigate to dataset directory
cd workflow/datasets/novogene_llc/

# Dry-run (preview DAG)
snakemake --cores 8 -np

# Run full pipeline
snakemake --cores 8

# Run single phase
snakemake phase08_done --cores 8

# Visualize DAG
snakemake --dag | dot -Tpdf > pipeline_dag.pdf
```

### **B. Manual Execution (Development/Debugging)**
```bash
# Set environment variables
export DS_CONFIG="configs_private/datasets/INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1.yaml"
export DS_PREFIX="novogene_llc"

# Run phase orchestrator
Rscript workflow/scrna/08_de/run_de_pipeline.R

# Or run individual script
Rscript workflow/scrna/08_de/01_pseudobulk_de.R
```

### **C. CNV Module (Isolated Environment)**
```bash
# Switch to CNV module
cd modules/cnv/

# Set project root
export LUNGMETA_ROOT="/data1/morgen/biohub/projects/Lung_Cancer_Meta_2026_morgen"

# Run SCEVAN
Rscript scripts/01_scevan.R

# Return to main project
cd ../..

# Integrate CNV results
export DS_PREFIX="novogene_llc"
Rscript workflow/scrna/07_cnv/01_cnv_to_annotation.R
```

---

## **4. Dataset Config Conventions**

### **Config Priority Chain (3 Layers)**
```
Dataset Private Config (configs_private/datasets/{id}.yaml)  ← Highest priority
    ↓ overrides
Project Params (configs/params/scrna_qc_params.yaml)         ← Mid priority
    ↓ overrides
Hardcoded Defaults (in R functions)                          ← Fallback
```

### **Dataset Config Structure**
```yaml
# configs_private/datasets/INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1.yaml
dataset_id: "INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1"
species: "mouse"
n_samples: 6

groups:
  mc:  { description: "empty vector control", samples: ["mc_1", "mc_2"] }
  FL:  { description: "wildtype overexpression", samples: ["FL_1", "FL_2"] }
  A1:  { description: "mutant overexpression", samples: ["A1_1", "A1_2"] }

samples:
  A1_1: { group: "A1", replicate: 1, fastq_prefixes: ["A1_1-1", "A1_1-2"] }
  # ... (other samples)

cellranger_out: "/data1/morgen/data_center/Lung_Cancer_2026/raw/internally_generated/INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1/cellranger_out"

# Optional: dataset-specific QC overrides
qc_overrides:
  n_mad: 5  # Override project default (3)
```

### **Species-specific Configs**
- **Mouse**: `configs/params/scrna_annotation_params.yaml`
  - SingleR refs: `celldex::MouseRNAseqData`, `celldex::ImmGenData`
  - Mito pattern: `^mt-`
- **Human**: `configs/params/scrna_annotation_params_human.yaml`
  - SingleR refs: `celldex::HumanPrimaryCellAtlasData`, `celldex::BlueprintEncodeData`
  - Mito pattern: `^MT-`

### **DE Contrasts (Project-specific)**
- **INT_Novogene**: `configs/annotation/de_contrasts.yaml`
  - 3 contrasts: FL_vs_mc, A1_vs_mc, A1_vs_FL
  - 3 axes: tumor_intrinsic, tme_per_celltype, epithelial_normal
- **GSE253718**: `configs/annotation/de_contrasts_gse253718.yaml`
  - 1 contrast: TKI_vs_naive
  - 2 axes: tumor_intrinsic, tme_per_celltype

---

## **5. scRNA Phase Definitions**

| **Phase** | **Name** | **Input** | **Output** | **Status** |
|-----------|----------|-----------|-----------|-----------|
| **01** | Alignment | FASTQ | `filtered_feature_bc_matrix/` | ✅ Done |
| **02** | QC | Cell Ranger output | `clean/{sid}_clean.rds` | ✅ Done |
| **03** | Normalize | Per-sample clean RDS | `pca/merged_pca.rds` | ✅ Done |
| **04** | Integrate | PCA object | `harmony/seurat_harmony.rds` | ✅ Done |
| **05** | Cluster | Harmony object | `objects/seurat_clustered.rds` | ✅ Done |
| **06** | Annotate | Clustered object + CSV | `objects/seurat_annotated_final.rds` ★★ | ✅ Done |
| **07** | CNV | Annotated object | CNV labels in metadata | ✅ Done |
| **08** | DE | Final annotated object | `deseq2/*.csv`, `enrichment/*.csv` | ✅ Done |
| **09** | CellChat | Final annotated object | Cell-cell communication networks | 🔴 Next |
| **10** | Trajectory | Final annotated object | Lineage inference | ⬜ Planned |
| **11** | Advanced | Final annotated object | SCENIC/Velocity/Scoring | ⬜ Planned |

### **Phase Details**

See `workflow/scrna/README.md` for detailed phase descriptions.

**Key Object**: `results/scrna/{ds_prefix}/06_annotate/objects/seurat_annotated_final.rds` is the **Single Source of Truth** for all downstream analyses (Phase 08-11).

---

## **6. Common Commands and Expected Outputs**

### **Check Environment Status**
```bash
# R environment
Rscript -e 'renv::status()'

# Snakemake environment
conda activate snakemake_env
snakemake --version
```

### **Run Full Pipeline (Snakemake)**
```bash
cd workflow/datasets/novogene_llc/
snakemake --cores 16 2>&1 | tee ../../logs/snakemake_$(date +%Y%m%d_%H%M%S).log
```
**Expected outputs**:
- Sentinel files: `results/scrna/novogene_llc/{phase}/.phase{N}_done`
- Final object: `results/scrna/novogene_llc/06_annotate/objects/seurat_annotated_final.rds` (13.5 GB)
- DE results: `results/scrna/novogene_llc/08_de/deseq2/all_de_results_combined.csv` (125 MB)

### **Approve Quality Gates**
```bash
# After reviewing Phase 05 clustering
touch results/scrna/novogene_llc/05_cluster/.gate_approved

# After reviewing Phase 06 annotation
touch results/scrna/novogene_llc/06_annotate/.gate_approved
```

### **Debug Configuration**
```bash
# Check final QC params
Rscript -e 'source("workflow/scrna/functions/qc_utils.R"); print(load_qc_params())'

# Check dataset config
Rscript -e 'source("workflow/scrna/functions/io_scrna.R"); Sys.setenv(DS_CONFIG="configs_private/datasets/INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1.yaml"); print(load_dataset_config())'
```

---

## **7. Rules for Adding a New Dataset**

### **Quick Steps**
1. Prepare data on HPC
2. Create public config: `configs/datasets/{dataset_id}.yaml`
3. Create private config: `configs_private/datasets/{dataset_id}.yaml`
4. Create DE contrasts: `configs/annotation/de_contrasts_{project}.yaml` (if needed)
5. Create Snakemake entry: `workflow/datasets/{short_name}/`
6. Register dataset: Add to `metadata/registry/dataset_registry.csv`
7. Run pipeline: `cd workflow/datasets/{short_name}/ && snakemake --cores 16`

**Detailed instructions**: See `docs/SOP_DATASET_ONBOARDING.md`

---

## **8. Rules for Debugging Pipeline Failures**

### **Common Failure Modes**

#### **1. Missing Input File**
```
Error: Input not found: results/scrna/novogene_llc/05_cluster/objects/seurat_clustered.rds
```
**Solution**: Check if previous phase completed, re-run if needed.

#### **2. Environment Variable Not Set**
```
Error: DS_PREFIX not set
```
**Solution**:
```bash
export DS_CONFIG="configs_private/datasets/{dataset_id}.yaml"
export DS_PREFIX="{short_name}"
```

#### **3. Config Priority Confusion**
**Solution**: Debug with `load_qc_params()` to see final effective config.

#### **4. CNV Output Contract Violation**
```
Error: undefined columns: scevan_label
```
**Solution**: Check CNV module output format, verify required columns exist.

#### **5. Celltype Mapping Out of Sync**
```
⚠️  Unmapped clusters: 27, 28, 29, 30
```
**Solution**: Update `configs/annotation/celltype_mapping.csv` with new clusters.

**More details**: See Section 8 in full document.

---

## **9. Naming Conventions**

### **Dataset IDs**
```
{SOURCE}_{VENDOR/ACCESSION}_{DISEASE}_{MODALITY}_{DESCRIPTION}_{DATE}_{VERSION}

Examples:
- INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1
- PUB_GEO_LUAD_scRNA_GSE253718_v1
```

### **Dataset Short Names (ds_prefix)**
```
{project}_{species/model}

Examples:
- novogene_llc
- gse253718
```

### **Metadata Columns**
- **Core**: `sample_id`, `group`, `nFeature_RNA`, `nCount_RNA`, `percent.mt`
- **Cell cycle**: `S.Score`, `G2M.Score`, `Phase`, `CC.Difference`
- **Clustering**: `seurat_clusters`, `clusters_res0.3`, etc.
- **Annotation**: `celltype_L1`, `celltype_L2`, `annotation_confidence`, `is_artifact`
- **CNV**: `scevan_label`, `scevan_subclone`

---

## **10. Known Ambiguities / TODOs**

### **Need Confirmation**
1. Phase 09-11 priority and target cell types
2. Meta-analysis strategy (cross-dataset integration method)
3. Template-ization timeline
4. CNV module version control mechanism

### **Known Issues**
1. Backup files (`.bak_20260409`) should be removed
2. Documentation gaps: CNV contract, config priority guide, function library docs
3. Test datasets in registry (9 `VEND_Test_*` entries)
4. Large object fragility (13.5 GB RDS, no backup/checksum)

### **Future Enhancements**
1. Contract validation functions
2. Environment check in all orchestrators
3. Quick runner script for manual execution
4. Automated testing for config chain and CNV contract

---

## **11. New Collaborator Onboarding: Shortest Path**

### **Day 1: Understand the Big Picture (30 min)**
1. Read `README.md`
2. Read this document (sections 1-5)
3. Explore `results/scrna/novogene_llc/`

### **Day 2: Setup Environment (1-2 hours)**
1. Read `docs/SETUP.md`
2. Create data symlink: `ln -s /data1/morgen/data_center/Lung_Cancer_2026 data`
3. Restore R environment: `Rscript -e 'renv::restore()'`
4. Activate Snakemake: `conda activate snakemake_env`

### **Day 3: Run a Test Phase (2-3 hours)**
1. Re-run Phase 08 DE manually
2. Compare outputs with existing results
3. Verify plots

### **Day 4: Understand Config System (1-2 hours)**
1. Read Section 4 (Dataset Config Conventions)
2. Debug config priority with `load_qc_params()`
3. Compare project vs dataset configs

### **Day 5: Add a New Dataset (Practice)**
1. Read Section 7 (Rules for Adding a New Dataset)
2. Practice with GSE253718 (already exists)
3. Run Phase 02 QC on one sample

---

## **12. Quick Reference: File Paths**

### **Critical Objects**
```
# Single Source of Truth
results/scrna/{ds_prefix}/06_annotate/objects/seurat_annotated_final.rds

# DE results
results/scrna/{ds_prefix}/08_de/deseq2/all_de_results_combined.csv
results/scrna/{ds_prefix}/08_de/enrichment/all_enrichment_combined.csv
```

### **Key Config Files**
```
# Dataset configs
configs/datasets/{dataset_id}.yaml                           # Public
configs_private/datasets/{dataset_id}.yaml                   # Private

# Analysis params
configs/params/scrna_qc_params.yaml                          # QC
configs/params/scrna_annotation_params.yaml                  # Mouse annotation
configs/params/scrna_annotation_params_human.yaml            # Human annotation

# Manual annotation
configs/annotation/celltype_mapping.csv

# DE contrasts
configs/annotation/de_contrasts.yaml                         # INT_Novogene
configs/annotation/de_contrasts_gse253718.yaml               # GSE253718
```

### **Utility Functions**
```
# Layer 1: Project-wide
scripts/utils/utils_io.R
scripts/utils/utils_plotting.R
scripts/utils/utils_registry.R

# Layer 2: scRNA-specific
workflow/scrna/functions/io_scrna.R
workflow/scrna/functions/qc_utils.R
workflow/scrna/functions/annotation_utils.R
```

---

## **13. Contact & Collaboration**

### **Primary Analyst**
- **Name**: morgen
- **Role**: Pipeline architect, primary analyst
- **Focus**: scRNA-seq pipeline development, INT_Novogene analysis

### **Collaboration Protocol**
1. Before modifying code: Discuss design
2. Before adding dataset: Follow Section 7
3. Before changing config: Understand priority chain (Section 4)
4. After major changes: Update this document and README.md

---

## **14. Common Pitfalls**

### **Pitfall 1: Running Scripts Without Environment Variables**
```bash
# ❌ Wrong
Rscript workflow/scrna/08_de/run_de_pipeline.R

# ✅ Correct
export DS_CONFIG="configs_private/datasets/INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1.yaml"
export DS_PREFIX="novogene_llc"
Rscript workflow/scrna/08_de/run_de_pipeline.R
```

### **Pitfall 2: Config Override Not Recognized**
Dataset config overrides project config. Check `qc_overrides` in dataset YAML.

### **Pitfall 3: Forgetting Quality Gates**
Snakemake stops at gates. Review and approve: `touch .gate_approved`

### **Pitfall 4: Celltype Mapping Out of Sync**
After re-clustering, update `celltype_mapping.csv` with new clusters.

### **Pitfall 5: CNV Module Format Changed**
If CNV output columns change, update `07_cnv/01_cnv_to_annotation.R`.

---

## **15. Glossary**

| **Term** | **Definition** |
|----------|---------------|
| **ds_prefix** | Short dataset name (e.g., `novogene_llc`) |
| **dataset_id** | Full identifier (e.g., `INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1`) |
| **Phase** | Major analysis stage (01-11) |
| **Sentinel file** | `.phase{N}_done` marker |
| **Quality gate** | Manual review checkpoint (`.gate_approved`) |
| **Layer 1/2/3** | Function hierarchy (project/scRNA/phase-specific) |
| **MAD** | Median Absolute Deviation (adaptive QC) |
| **Pseudobulk** | Sample-level aggregation for DE |
| **celltype_L1/L2** | Cell type hierarchy (broad/specific) |
| **CNV** | Copy Number Variation |
| **Contract** | Interface specification (input/output format) |

---

## **16. Emergency Contacts & Resources**

### **When Things Go Wrong**

| **Issue** | **First Action** | **Escalation** |
|-----------|-----------------|----------------|
| Pipeline fails | Check Section 8 | Review logs in `logs/snakemake/` |
| Config not working | Debug with `load_qc_params()` | Check Section 4 |
| Object corrupted | Check file size, try `readRDS()` | Restore backup or re-run |
| Out of disk space | Check `du -sh results/` | Clean or expand storage |
| R package missing | `renv::restore()` | Check `renv.lock` |
| Snakemake DAG error | `snakemake -np` | Check Snakefile syntax |

### **Key Documentation**
- `README.md` — Project overview
- `docs/SETUP.md` — Environment setup
- `docs/SOP_DATASET_ONBOARDING.md` — Dataset onboarding
- `workflow/scrna/README.md` — Pipeline architecture
- This document — Working memory

---

## **17. Version History**

| **Version** | **Date** | **Changes** | **Author** |
|-------------|----------|-------------|-----------|
| 1.0.0 | 2026-04-20 | Initial Working Memory created | AI Assistant + morgen |

---

**END OF PROJECT WORKING MEMORY v1.0.0**

---

## **Related Documents**

- **Setup & Environment**: `docs/SETUP.md`
- **Dataset Onboarding**: `docs/SOP_DATASET_ONBOARDING.md`
- **Pipeline Architecture**: `workflow/scrna/README.md`
- **Project Overview**: `README.md`
