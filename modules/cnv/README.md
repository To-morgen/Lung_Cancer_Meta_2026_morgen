# Module: CNV Analysis (SCEVAN / inferCNV)

**结论先行 (BLUF)**
本模块是单细胞转录组流水线中的**独立计算节点 (Isolated Compute Node)**。其核心设计初衷是通过 **环境隔离 (Environment Isolation)** 技术，解决 `SCEVAN` 与 `inferCNV` 等重型生信工具在 R 语言环境下的依赖冲突 (Dependency Conflicts)，确保主流水线的简洁性与稳健性。

---

## Purpose (核心目标)

* **隔离性 (Isolation):** 专用于从 scRNA-seq 数据中推断拷贝数变异 (CNV)。
* **兼容性保障:** 鉴于 `SCEVAN` 等工具依赖链条极其复杂且脆弱 (Heavy/Fragile Dependencies)，将其从主项目 `renv` 环境中剥离，防止污染全局依赖树。

---

## Environment (运行环境)

该模块在 `modules/cnv/` 目录下维护一套完全独立的 R 包生态系统。



| 关键组件 | 说明 |
| :--- | :--- |
| **Package Manager** | `renv` (独立锁文件 `renv.lock`) |
| **Core Algorithms** | `SCEVAN`, `inferCNV` |
| **Viz & Math** | `ComplexHeatmap`, `circlize`, `rjags` (JAGS 贝叶斯推断) |
| **Dependency Risk** | 高 (依赖系统级库文件，建议在固定镜像或独立 Conda 环境下执行) |

---

## Usage (执行指南)

在执行 CNV 计算时，需手动切换上下文至模块目录，以激活正确的 `renv` 环境。

```bash
# 进入 CNV 专用模块目录
cd modules/cnv

# 1. 恢复独立的 R 环境
Rscript -e 'renv::restore()'

# 2. 执行核心算法 (SCEVAN Inference)
Rscript scripts/01_scevan.R

# 3. 执行结果验证与 QC
Rscript scripts/02_validate_scevan.R
```

---

## Integration (跨目录整合)

CNV 分析完成后，计算出的恶性标签 (Malignant Labels) 需要被“拉回”到主分析流水线中进行细胞类型重注释。

```bash
# 返回项目根目录 (Main Project Root)
cd ~/biohub/projects/Lung_Cancer_Meta_2026_morgen

# 执行回流脚本：将模块输出整合至主项目 Seurat 对象
Rscript workflow/scrna/07_cnv/01_cnv_to_annotation.R
```

---
