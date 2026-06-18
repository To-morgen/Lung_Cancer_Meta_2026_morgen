# Lung Cancer Meta Pipeline

A reproducible Snakemake + R pipeline for lung cancer single-cell meta-analysis.

This repository is intended to be a public, reusable pipeline asset: it stores code,
sanitized configuration templates, and documentation. Raw data, large intermediate
objects, generated results, private handoff notes, and local agent/tooling state are
intentionally excluded.

## What the pipeline covers

```text
Dataset config
  -> 01 alignment / count-matrix intake
  -> 02 QC
  -> 03 normalization
  -> 04 integration
  -> 05 clustering
  -> 06 annotation
  -> 07 CNV review
  -> 08 pseudobulk DE / enrichment
  -> 09 CellChat communication analysis
```

The current public milestone is **v0.1**: a public-safe, config-driven scaffold for
lung cancer scRNA dataset onboarding and downstream analysis.

## Repository scope

Tracked in Git:

- reusable R scripts and Snakemake rules
- public-safe YAML/CSV templates
- module environment lockfiles
- setup and onboarding documentation

Not tracked in Git:

- raw FASTQ/BAM/count matrices
- Seurat, AnnData, inferCNV, CellChat, or other large objects
- generated plots, reports, logs, and Snakemake runtime state
- private dataset configs, local notes, handoff files, and machine-specific paths

## Layout

```text
configs/
  annotation/          # annotation maps, CNV target templates, DE contrast templates
  datasets/            # dataset metadata templates and public dataset configs
  params/              # shared analysis parameters
modules/
  cnv/                 # SCEVAN / inferCNV scripts in an isolated renv
  cellchat/            # CellChat scripts in an isolated renv
workflow/
  datasets/            # per-dataset Snakemake entry points
  scrna/               # phase-level R scripts
  snakemake/rules/     # reusable Snakemake rules
docs/                  # public setup and onboarding docs
```

## Quick start

```bash
# Clone
git clone <repo_url>
cd Lung_Cancer_Meta_2026_morgen

# Copy and edit local environment hints
cp .env.sh.example .env.sh

# Restore the main R environment
Rscript -e 'renv::restore()'

# Restore module environments when needed
(cd modules/cnv && Rscript -e 'renv::restore()')
(cd modules/cellchat && Rscript -e 'renv::restore()')
```

## Run a dataset workflow

Each dataset workflow has its own `workflow/datasets/<dataset>/Snakefile` and
`config.yaml`. Start with a dry-run.

```bash
snakemake -s workflow/datasets/example_lung_2grp/Snakefile -n --cores 1
```

For a new private dataset, create local files from the templates and keep them out
of Git:

```bash
cp configs/datasets/_TEMPLATE_internally_generated_.yaml \
  configs_private/datasets/<dataset_id>.yaml

cp workflow/datasets/example_lung_2grp/config.yaml \
  workflow/datasets/<your_dataset>/config.yaml
```

Use relative paths through `data/` or placeholders in public configs. Put real raw
data paths, vendor delivery paths, and machine-specific settings in `configs_private/`
or `.env.sh`.

## Public v0.1 focus

Milestone v0.1 focuses on the public pipeline surface:

- dataset-driven Snakemake orchestration
- manual gate rules for review checkpoints
- CNV module integration with SCEVAN/inferCNV scoring
- pseudobulk DE configuration
- CellChat prepare/run/compare/deep-dive/summary scripts
- public-safe ignore rules and setup docs

## Design principles

1. **Config-driven analysis** — dataset, contrast, threshold, and path choices live in YAML/CSV files.
2. **Code-data separation** — Git stores pipeline logic; data and generated outputs live outside Git.
3. **Module isolation** — dependency-heavy tools such as CNV and CellChat use module-local renv environments.
4. **Manual gates for biological review** — annotation and downstream interpretation require explicit review checkpoints.
5. **Public/private separation** — public code is reusable; private analysis notes and unpublished results stay local/private.

## Documentation

- `docs/pipeline_overview.md` — phase map, repository contracts, and v0.1 scope
- `docs/dataset_onboarding.md` — public-safe process for adding new datasets
- `docs/SETUP.md` — setup, storage policy, and execution conventions
- `docs/SOP_DATASET_ONBOARDING.md` — dataset onboarding lifecycle
- `configs/datasets/_template.yaml` — public dataset metadata template
- `configs/datasets/_TEMPLATE_internally_generated_.yaml` — private dataset config template

## License and citation

This repository is under active development. Add a project license and citation
file before relying on it as a stable external dependency.
