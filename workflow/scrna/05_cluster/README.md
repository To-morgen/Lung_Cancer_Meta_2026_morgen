# Phase 5: Clustering + UMAP

## Scripts
1. `01_cluster_umap.R` — FindNeighbors + FindClusters (multi-res) + RunUMAP
2. `02_cluster_qc.R` — Cluster composition, QC metrics, clustree

## Key Parameters
- Reduction: `harmony` (NOT raw PCA)
- Algorithm: Leiden (4)
- Resolutions tested: 0.3, 0.5, 0.8, 1.0, 1.2
- Default: 0.8

## Decision Point ✋
After running, review:
1. `04_clustree.png` — resolution stability
2. `05_cluster_composition.png` — balanced sample contribution?
3. `02_multi_resolution.png` — which resolution gives biologically meaningful clusters?

Then update `scrna_qc_params.yaml → clustering.default_resolution` if needed.

## Output
results/scrna/05_cluster/
├── objects/ seurat_clustered.rds (all resolutions stored)
├── plots/ 01-07 diagnostic plots
└── reports/ composition tables, resolution comparison

