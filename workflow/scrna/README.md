# scRNA-seq Analysis Workflow

## Complete Pipeline

阶段一: 单样本清洗 (per-sample)
────────────────────────────────────────
01_alignment/ Cell Ranger count ✅ Done
02_qc/ SoupX → Seurat → MAD QC 🔴 Current
→ scDblFinder → clean

阶段二: 汇集与特征提取
────────────────────────────────────────
03_normalize/ merge → SCTransform ⬜ Next
→ Cell Cycle Scoring
→ PCA

阶段三: 去批次与聚类
────────────────────────────────────────
04_integrate/ Harmony ⬜ Planned
05_cluster/ FindNeighbors/FindClusters ⬜ Planned
→ UMAP

阶段四: 生物学注释
────────────────────────────────────────
06_annotate/ SingleR (auto) ⬜ Planned
→ FindAllMarkers
→ Manual annotation

06_downstream/ DE / Pathway / CellChat ⬜ Future


## Key Design Decisions

1. **Doublet removal BEFORE merge** — physical artifact, must be per-sample
2. **MAD-based QC** — data-driven, n_mad=3, log space for count metrics
3. **SCTransform** over LogNormalize — better variance stabilization
4. **Cell Cycle**: Score always, regress only if needed (CC.Difference strategy for tumor study)
5. **Harmony** for batch correction — fast, effective, works in PCA space
6. **Leiden** clustering (algorithm 4) with multiple resolutions

## Cell Cycle Strategy (for this tumor project)

Always score S.Score + G2M.Score (CellCycleScoring)
Run SCTransform WITHOUT regression first
Check PCA/UMAP:
If G2M/S dominate PC1/PC2 → regress CC.Difference
If minimal effect → keep as-is
For LLC lung cancer: proliferation is key biology
→ Prefer CC.Difference (preserves cycling vs quiescent)
→ Never blindly remove S.Score + G2M.Score

## Function Layers

Layer 1: scripts/utils/ Project-wide
Layer 2: workflow/scrna/functions/ scRNA-specific
Layer 3: workflow/scrna/02_qc/*.R Step scripts

