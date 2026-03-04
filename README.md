# Lung Cancer Multi-omics Atlas 2026
Publication-grade multi-omics integration for lung cancer, focusing on **curation**, **reproducible pipelines**, and **analysis-ready outputs**.

## Project structure
- `data/` : symlink to local data center (NOT tracked as files). See `docs/SETUP.md`.
- `metadata/` : sample sheets, manifests, and study-level annotations
- `workflow/` : Snakemake pipelines
- `scripts/` : analysis scripts (R/Python)
- `configs/` : configs for pipelines and tools
- `results/` : derived results (ignored by default; keep lightweight summaries if needed)
- `figures/` : final figures (ignored by default; keep exported panels separately if desired)
- `docs/` : method cards / notes
- `logs/`, `tmp/` : runtime artifacts (ignored)

## Reproducibility
Environment files will be provided:
- `environment.yml` (conda/mamba) and/or `renv.lock` (R)

## Data governance
Raw patient-level data are not uploaded. This repository tracks:
- acquisition scripts, processing code, metadata, and manifests.

## Setup
See `docs/SETUP.md` for HPC setup (data symlink, governance rules).
