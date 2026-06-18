# Phase 07: CNV Integration

**结论先行 (BLUF)**
Phase 07 是整个单细胞流水线的“恶性鉴定 (Malignancy Proof)”阶段。通过在独立环境 (Isolated Environment) 中运行 **SCEVAN/inferCNV**，将染色质拷贝数变异 (CNV) 证据回传至主分析对象，从而实现从“表型聚类 (Phenotypic Clustering)”到“基因型验证 (Genotypic Validation)”的跨越。

---

## Architecture (架构设计)

为了解决生信分析中常见的“依赖冲突 (Dependency Hell)”，CNV 分析被物理拆分为两个执行环境：

| Component | Function | Location | R Environment (renv) |
| :--- | :--- | :--- | :--- |
| **Computation** | SCEVAN/inferCNV 核心算法运行 | `modules/cnv/scripts/` | **modules/cnv** (独立) |
| **Integration** | 标签回流 (Label Backfilling) 与验证 | `workflow/scrna/07_cnv/` | **main project** (主环境) |

---

## Why split? (环境拆分逻辑)

* **依赖解耦 (Dependency Decoupling):** SCEVAN 与 inferCNV 引入了极其沉重的底层依赖（如 `rjags`, `ComplexHeatmap`, `HiddenMarkov`）。这些包的版本要求往往与主分析环境中的 `Seurat v5` 或 `SCTransform` 存在冲突。
* **环境轻量化:** 避免主环境 `renv.lock` 过于庞大，提升服务器上的环境构建速度。

---

## Data Flow (数据流向)



```text
1. 核心计算 (Module)
   modules/cnv/scripts/01_scevan.R
     └── Output: results/scrna/07_cnv/scevan/seurat_with_scevan.rds

2. 标签整合 (Main Workflow)
   workflow/scrna/07_cnv/01_cnv_to_annotation.R
     ├── Reads: seurat_with_scevan.rds (来自计算层)
     ├── Reads: celltype_mapping.csv (来自配置层)
     └── Writes: seurat_annotated_final.rds (最终对象)
```

---

## Execution Workflow (执行流程)

### Step 1: Module Environment (计算层)
在该步骤中，利用 `modules/cnv` 目录下特定的 `renv` 执行高负荷计算。

```bash
# 进入模块目录
cd modules/cnv

# 运行 SCEVAN 算法
Rscript scripts/01_scevan.R

# 运行验证逻辑 (QC for CNV)
Rscript scripts/02_validate_scevan.R
```

### Step 2: Main Project Integration (整合层)
回到主项目环境，将计算得到的 CNV 结果（如 `Malignant` vs `Normal` 标签）合并到最终的 Seurat 对象中。

```bash
# Return to the main project root
cd /path/to/Lung_Cancer_Meta_2026_morgen

# 执行标签回流脚本
Rscript workflow/scrna/07_cnv/01_cnv_to_annotation.R
```

---
