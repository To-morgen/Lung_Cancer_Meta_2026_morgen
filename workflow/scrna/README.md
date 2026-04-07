# scRNA-seq Analysis Workflow

**Project Context:** Multi-omics Atlas 2026  
**Primary Runtime:** R / Seurat Ecosystem (Ubuntu HPC)  
**Workflow Engine:** Config-driven (YAML/CSV)

-----

## Pipeline Architecture

| Phase | Description | Environment | Status |
| :--- | :--- | :--- | :--- |
| **Phase 1: Per-sample** | `01_alignment`: Cell Ranger Count<br>`02_qc`: SoupX → Seurat → MAD QC → scDblFinder | main renv | ✅ Done |
| **Phase 2: Integration** | `03_normalize`: Merge → CC Score → SCTransform → PCA<br>`04_integrate`: Harmony batch correction<br>`05_cluster`: Leiden clustering → UMAP | main renv | ✅ Done |
| **Phase 3: Annotation** | `06_annotate`: FindMarkers → SingleR → Manual<br>`07_cnv`: **SCEVAN** → Validate → **回流** (Refine labels) | mixed renv | ✅ Done |
| **Phase 4: Downstream** | `08_de`: Pseudobulk DE (DESeq2) + GSEA<br>`09_cellchat`: Cell communication (L1 first)<br>`10_trajectory`: Monocle3 / Slingshot | main renv | 🔴 Next |
| **Phase 5: Advanced** | `11_advanced`: SCENIC / Velocity / Scoring | modules/ | ⬜ Planned |

-----

## Object Version Chain

The following flow maintains the integrity of the data lineage, ensuring that CNV-validated tumor labels are propagated to all downstream analyses.

1.  **Input Cluster Object** `05_cluster/objects/seurat_clustered.rds`
2.  **Initial Annotation (v1: pre-CNV)** `06_annotate/objects/seurat_annotated.rds`  
    *Primary lineage assignment via markers/SingleR.*
3.  **Final Annotated Object (v2: post-CNV)** `06_annotate/objects/seurat_annotated_final.rds`  
    *Incorporates SCEVAN labels + Tumor\_putative confirmation.*
      * **Downstream Entry Point:** This object is the **Single Source of Truth** for:
          * `08_de/` (Differential Expression)
          * `09_cellchat/` (Interactome)
          * `10_trajectory/` (Lineage Inference)

-----

## Analysis Design

### Axis 1: Tumor-intrinsic transcriptional differences

  * **Comparison:** FL vs A1 vs mock (within confirmed **Tumor** clusters).
  * **Method:** Pseudobulk **DESeq2** → **fgsea** (Hallmark, KEGG, Reactome, GO\_BP).

### Axis 2: TME composition & communication differences

  * **Comparison:** FL vs A1 vs mock (per immune/stromal **celltype\_L2**).
  * **Method:** Pseudobulk DE + Compositional analysis (**sccomp**) + **CellChat**.

-----

## Configuration Management

  * **`configs/params/scrna_qc_params.yaml`**: Defines MAD thresholds ($n\_mad$) and hard cutoffs.
  * **`configs/params/scrna_annotation_params.yaml`**: Stores lineage markers and project-specific gene sets.
  * **`configs/annotation/celltype_mapping.csv`**: Human-curated mapping of clusters to cell identities.
  * **`configs/annotation/de_contrasts.yaml`**: Specifies DE groups, methods, and enrichment targets.

-----

## Functions (Layer Architecture)

To ensure code maintainability, functions are divided into project-wide utilities and workflow-specific logic:

  * **Layer 1: Project-wide Utilities**
      * `scripts/utils/utils_io.R`: Generic Read/Write operations.
      * `scripts/utils/utils_plotting.R`: Standardized ggplot2 themes and palettes.
  * **Layer 2: Domain-specific Logic**
      * `workflow/scrna/functions/io_scrna.R`: Custom Seurat object loaders and path management.
      * `workflow/scrna/functions/qc_utils.R`: MAD calculation and metadata cleaning.
      * `workflow/scrna/functions/annotation_utils.R`: Config-to-Seurat mapping and marker extraction.

-----
