# Workflow layer

This directory stores **project-wide operational workflows** for the Lung Cancer Multi-omics Atlas 2026 project.

It is designed for **data governance and standardization**, not for downstream scientific analysis modules.

## Structure

- `intake/`  
  Dataset registration, file manifest generation, and structural intake checks.

- `qc/`  
  Quality control workflows for different dataset classes.

- `harmonize/`  
  Standardization and alignment of data objects into internally curated formats.

- `promote/`  
  Promotion of datasets from curated objects to atlas-ready objects.

## Design principles

1. **General utilities** belong to `/scripts`
2. **Project-wide workflows** belong to `/workflow`
3. **Scientific analysis modules** belong to `/modules`
4. **Dataset-specific differences** should be expressed through `/configs/datasets/*.yaml`
5. Heavy data objects remain in private HPC-backed `data/`

## Current flow

Typical vendor-curated dataset flow:

`register -> manifest_check -> intake_qc -> harmonize -> to_curated -> atlas_ready`

Typical public-original dataset flow:

`register -> manifest_check -> raw_qc -> harmonize -> to_curated -> atlas_ready`