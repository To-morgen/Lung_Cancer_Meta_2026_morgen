# Phase 6: Cell Type Annotation

## Workflow

Step 1 (auto): FindAllMarkers → marker tables + heatmap + lineage dotplot + GPRIN2 plot
Step 2 (auto): SingleR → MouseRNAseqData + ImmGenData → per-cluster + per-cell labels
─── HUMAN REVIEWS RESULTS AND FILLS CSV ───
Step 3 (manual): Apply CSV mapping → seurat_annotated.rds (canonical output)


## How to use

```bash
# 1. Run automated steps
source .env.sh
Rscript workflow/scrna/06_annotate/run_annotation_pipeline.R 2>&1 | tee logs/phase6_auto.log

# 2. Review results
#    - plots/01-05: marker heatmap, dotplots, lineage markers, GPRIN2
#    - plots/06-08: SingleR UMAP + score heatmaps
#    - markers/top10_markers_per_cluster.csv
#    - reports/singler_cluster_annotation.csv

# 3. Create annotation mapping
cp configs/annotation/celltype_mapping_template.csv configs/annotation/celltype_mapping.csv
# Edit celltype_mapping.csv with your annotations

# 4. Apply manual annotation
Rscript workflow/scrna/06_annotate/run_annotation_pipeline.R --manual 2>&1 | tee logs/phase6_manual.log
Annotation Strategy: Two-tier (L1 + L2)
L1 (coarse):  Tumor / T_cell / Myeloid / B_cell / NK / Stroma / Cycling / Unknown
L2 (fine):    LLC_main / CD8_Tex / TAM_M2 / DC_cDC1 / ...
Key Output
results/scrna/06_annotate/
├── objects/
│   ├── seurat_singler_annotated.rds    (intermediate: with SingleR labels)
│   └── seurat_annotated.rds            (CANONICAL: all downstream reads this)
├── markers/
│   ├── all_markers_res0.8.csv
│   ├── top{3,5,10,20}_markers_per_cluster.csv
├── singler/
│   └── singler_*.rds
├── plots/
│   ├── 01-05: marker/lineage/GPRIN2 plots
│   ├── 06-08: SingleR diagnostic plots
│   └── 09-11: final annotated UMAP + composition
└── reports/
    ├── singler_cluster_annotation.csv
    ├── composition_L1.csv
    ├── composition_L2.csv
    ├── composition_by_group.csv
    └── cell_metadata_annotated.csv
