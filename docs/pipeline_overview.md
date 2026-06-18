# Pipeline Overview

Lung Cancer Meta Pipeline is a config-driven Snakemake + R workflow for lung
cancer single-cell RNA-seq analysis. The public repository focuses on reusable
pipeline logic rather than private data or unpublished results.

## High-level flow

```text
Dataset config
  │
  ├─ 01 alignment / count-matrix intake
  │    FASTQ or existing 10x matrices -> per-sample count matrices
  │
  ├─ 02 QC
  │    SoupX / Seurat object creation / MAD filters / doublet removal
  │
  ├─ 03 normalization
  │    merge samples / cell-cycle scoring / SCTransform / PCA
  │
  ├─ 04 integration
  │    batch correction and integrated embeddings
  │
  ├─ 05 clustering
  │    UMAP / Leiden clustering / cluster QC
  │
  ├─ 06 annotation
  │    marker discovery / SingleR / reviewed cell-type mapping
  │
  ├─ 07 CNV
  │    SCEVAN / inferCNV scoring / annotation feedback
  │
  ├─ 08 DE
  │    pseudobulk DESeq2 / enrichment / plots
  │
  └─ 09 CellChat
       group-wise CellChat inference / comparison / summary reports
```

## Repository contracts

### Code-data separation

Git tracks pipeline logic and sanitized configuration templates. Data and outputs
are regenerated locally and stay outside Git.

```text
tracked:    workflow/, modules/, scripts/, configs/*templates*, docs/
ignored:    data/, results/, figures/, logs/, configs_private/, Notes/
```

### Dataset contract

Each dataset has a lightweight Snakemake entry point:

```text
workflow/datasets/<dataset>/
  Snakefile
  config.yaml
```

The Snakemake config declares the dataset ID, sample IDs, species, optional
phase-specific configs, and local input roots. Private paths belong in
`configs_private/` or `.env.sh`, not in public config files.

### Downstream object contract

The downstream phases expect an annotated Seurat object at:

```text
results/scrna/{ds_prefix}/06_annotate/objects/seurat_annotated_final.rds
```

This object is not tracked. It is the local handoff between annotation, CNV, DE,
and CellChat phases.

## Manual review gates

Some transitions are intentionally not fully automatic. For example, annotation
review and downstream approval require explicit gate files under `results/`.

This separates two different states:

- computational step completed
- biological interpretation approved

## Public v0.1 scope

v0.1 is intentionally narrow:

- public-safe README and setup docs
- example dataset workflow
- public dataset config example
- Phase01/07/08/09 Snakemake rules
- CNV module scripts and shared config
- CellChat summary/reporting script

Future milestones should add runnable toy data, automated smoke tests, and a
stable release tag after the public branch is reviewed.
