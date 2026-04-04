# Phase 4: Integration (Harmony)

## Script
- `01_harmony.R` — Harmony batch correction in PCA space

## Parameters (from scrna_qc_params.yaml)
- `integration.group_by`: `sample_id` (6 batches)
- `integration.n_pcs`: 30

## Output
results/scrna/04_integrate/
├── harmony/ seurat_harmony.rds
├── plots/ PCA pre/post comparison
└── reports/ batch centroid distances

