# **项目协作说明书：肺癌多组学图谱 2026**

**版本**: 1.0.0  
**最后更新**: 2026-04-20  
**分析师**: morgen  
**状态**: Phase 1-8 完成 (INT_Novogene), Phase 1-6 完成 (GSE253718)

---

## **1. 项目使命**

### **核心目标**
构建一个**配置驱动、可复现、模块化的 scRNA-seq 分析流程**，用于肺癌多组学研究，重点强调：
- 零硬编码参数（全部外置于 YAML/CSV）
- 基于阶段的架构（01-11）
- 多数据集支持（小鼠/人类，内部/公开）
- 依赖重的工具隔离环境（CNV, CellChat）
- 出版级可复现性（Level 3：外部研究者可复现）

### **当前焦点**
- **主数据集**: INT_Novogene_LLC（小鼠 LLC 模型，GPRIN2 野生型/突变体/对照，6 个样本）
- **验证数据集**: GSE253718（人类 LUAD，EGFR 突变，TKI 初治 vs 耐药，6 个样本）
- **分析阶段**: 完成 Phase 08（DE + 富集分析）后，进入 Phase 09（CellChat）

### **未来计划**
- 跨数据集的 meta 分析
- 流程模板化（流程成熟后）
- 全面的下游模块（CellChat, Trajectory, SCENIC, Velocity, Scoring）

---

## **2. 目录约定**

### **顶层结构**
```
Lung_Cancer_Meta_2026_morgen/
├── configs/                      # 公开配置（分析参数）
│   ├── annotation/               # 细胞类型映射、DE 对比设计
│   ├── datasets/                 # 公开数据集元数据（GEO 登录号）
│   ├── params/                   # QC、归一化、聚类参数
│   └── registry/                 # 数据集注册表 schema
├── configs_private/              # 私有配置（HPC 路径、未发表数据）
│   └── datasets/                 # 数据集特定路径和样本清单
├── data/                         # HPC 存储的符号链接（不追踪）
├── docs/                         # SOP、设置指南、流程文档
├── metadata/                     # 数据集注册表、接入报告
│   ├── registry/                 # dataset_registry.csv
│   └── reports/                  # 接入 QC 报告
├── modules/                      # 隔离的分析模块（独立 renv）
│   ├── cnv/                      # SCEVAN/inferCNV（独立 renv）
│   └── cellchat/                 # CellChat（计划中）
├── results/                      # 所有分析输出（不追踪）
│   └── scrna/{ds_prefix}/        # 每个数据集的结果
│       ├── 02_qc/ → 08_de/       # 阶段特定输出
│       └── .phase{N}_done        # 哨兵文件
├── scripts/                      # 项目级工具（Layer 1）
│   ├── setup/                    # 环境设置
│   └── utils/                    # IO、绘图、注册表辅助函数
├── workflow/                     # 分析流程
│   ├── datasets/{name}/          # 每个数据集的 Snakemake 入口
│   │   ├── Snakefile             # 数据集特定目标
│   │   └── config.yaml           # 样本、物种、路径
│   ├── scrna/                    # scRNA-seq 主流程
│   │   ├── 01_alignment/ → 11_advanced/  # 阶段目录
│   │   ├── functions/            # scRNA 专用工具（Layer 2）
│   │   └── {phase}/run_*_pipeline.R  # 阶段编排器
│   └── snakemake/rules/          # 共享 Snakemake 规则
├── renv.lock                     # 主 R 环境锁文件
└── README.md                     # 项目概览
```

### **关键原则**
- **Git 追踪**: 代码、配置、元数据、文档、renv.lock
- **Git 忽略**: data/, results/, logs/, .Rproj.user/, *.bak
- **HPC 存储**: 原始数据、处理后对象、分析输出
- **符号链接**: `data/` → `/data1/morgen/data_center/Lung_Cancer_2026/`

---

## **3. 主要工作流入口**

### **A. Snakemake（生产环境，推荐）**
```bash
# 导航到数据集目录
cd workflow/datasets/novogene_llc/

# 试运行（预览 DAG）
snakemake --cores 8 -np

# 运行完整流程
snakemake --cores 8

# 运行单个阶段
snakemake phase08_done --cores 8

# 可视化 DAG
snakemake --dag | dot -Tpdf > pipeline_dag.pdf
```

### **B. 手动执行（开发/调试）**
```bash
# 设置环境变量
export DS_CONFIG="configs_private/datasets/INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1.yaml"
export DS_PREFIX="novogene_llc"

# 运行阶段编排器
Rscript workflow/scrna/08_de/run_de_pipeline.R

# 或运行单个脚本
Rscript workflow/scrna/08_de/01_pseudobulk_de.R
```

### **C. CNV 模块（隔离环境）**
```bash
# 切换到 CNV 模块
cd modules/cnv/

# 设置项目根目录
export LUNGMETA_ROOT="/data1/morgen/biohub/projects/Lung_Cancer_Meta_2026_morgen"

# 运行 SCEVAN
Rscript scripts/01_scevan.R

# 返回主项目
cd ../..

# 整合 CNV 结果
export DS_PREFIX="novogene_llc"
Rscript workflow/scrna/07_cnv/01_cnv_to_annotation.R
```

---

## **4. 数据集配置约定**

### **配置优先级链（3 层）**
```
数据集私有配置 (configs_private/datasets/{id}.yaml)  ← 最高优先级
    ↓ 覆盖
项目参数 (configs/params/scrna_qc_params.yaml)         ← 中等优先级
    ↓ 覆盖
硬编码默认值 (在 R 函数中)                             ← 兜底
```

### **数据集配置结构**
```yaml
# configs_private/datasets/INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1.yaml
dataset_id: "INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1"
species: "mouse"
n_samples: 6

groups:
  mc:  { description: "空载对照", samples: ["mc_1", "mc_2"] }
  FL:  { description: "野生型过表达", samples: ["FL_1", "FL_2"] }
  A1:  { description: "突变体过表达", samples: ["A1_1", "A1_2"] }

samples:
  A1_1: { group: "A1", replicate: 1, fastq_prefixes: ["A1_1-1", "A1_1-2"] }
  # ... (其他样本)

cellranger_out: "/data1/morgen/data_center/Lung_Cancer_2026/raw/internally_generated/INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1/cellranger_out"

# 可选：数据集特定 QC 覆盖
qc_overrides:
  n_mad: 5  # 覆盖项目默认值（3）
```

### **物种特定配置**
- **小鼠**: `configs/params/scrna_annotation_params.yaml`
  - SingleR 参考: `celldex::MouseRNAseqData`, `celldex::ImmGenData`
  - 线粒体模式: `^mt-`
- **人类**: `configs/params/scrna_annotation_params_human.yaml`
  - SingleR 参考: `celldex::HumanPrimaryCellAtlasData`, `celldex::BlueprintEncodeData`
  - 线粒体模式: `^MT-`

### **DE 对比设计（项目特定）**
- **INT_Novogene**: `configs/annotation/de_contrasts.yaml`
  - 3 个对比: FL_vs_mc, A1_vs_mc, A1_vs_FL
  - 3 个分析轴: tumor_intrinsic, tme_per_celltype, epithelial_normal
- **GSE253718**: `configs/annotation/de_contrasts_gse253718.yaml`
  - 1 个对比: TKI_vs_naive
  - 2 个分析轴: tumor_intrinsic, tme_per_celltype

---

## **5. scRNA 阶段定义**

| **阶段** | **名称** | **输入** | **输出** | **状态** |
|-----------|----------|-----------|-----------|-----------|
| **01** | 比对 | FASTQ | `filtered_feature_bc_matrix/` | ✅ 完成 |
| **02** | QC | Cell Ranger 输出 | `clean/{sid}_clean.rds` | ✅ 完成 |
| **03** | 归一化 | 每样本 clean RDS | `pca/merged_pca.rds` | ✅ 完成 |
| **04** | 整合 | PCA 对象 | `harmony/seurat_harmony.rds` | ✅ 完成 |
| **05** | 聚类 | Harmony 对象 | `objects/seurat_clustered.rds` | ✅ 完成 |
| **06** | 注释 | 聚类对象 + CSV | `objects/seurat_annotated_final.rds` ★★ | ✅ 完成 |
| **07** | CNV | 注释对象 | 元数据中的 CNV 标签 | ✅ 完成 |
| **08** | DE | 最终注释对象 | `deseq2/*.csv`, `enrichment/*.csv` | ✅ 完成 |
| **09** | CellChat | 最终注释对象 | 细胞-细胞通讯网络 | 🔴 下一步 |
| **10** | 轨迹 | 最终注释对象 | 谱系推断 | ⬜ 计划中 |
| **11** | 高级 | 最终注释对象 | SCENIC/Velocity/Scoring | ⬜ 计划中 |

### **阶段详情**

详细阶段描述见 `workflow/scrna/README.md`。

**关键对象**: `results/scrna/{ds_prefix}/06_annotate/objects/seurat_annotated_final.rds` 是所有下游分析（Phase 08-11）的**唯一真值来源**。

---

## **6. 常用命令和预期输出**

### **检查环境状态**
```bash
# R 环境
Rscript -e 'renv::status()'

# Snakemake 环境
conda activate snakemake_env
snakemake --version
```

### **运行完整流程（Snakemake）**
```bash
cd workflow/datasets/novogene_llc/
snakemake --cores 16 2>&1 | tee ../../logs/snakemake_$(date +%Y%m%d_%H%M%S).log
```
**预期输出**:
- 哨兵文件: `results/scrna/novogene_llc/{phase}/.phase{N}_done`
- 最终对象: `results/scrna/novogene_llc/06_annotate/objects/seurat_annotated_final.rds` (13.5 GB)
- DE 结果: `results/scrna/novogene_llc/08_de/deseq2/all_de_results_combined.csv` (125 MB)

### **批准质量门控**
```bash
# 审查 Phase 05 聚类后
touch results/scrna/novogene_llc/05_cluster/.gate_approved

# 审查 Phase 06 注释后
touch results/scrna/novogene_llc/06_annotate/.gate_approved
```

### **调试配置**
```bash
# 检查最终 QC 参数
Rscript -e 'source("workflow/scrna/functions/qc_utils.R"); print(load_qc_params())'

# 检查数据集配置
Rscript -e 'source("workflow/scrna/functions/io_scrna.R"); Sys.setenv(DS_CONFIG="configs_private/datasets/INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1.yaml"); print(load_dataset_config())'
```

---

## **7. 添加新数据集的规则**

### **快速步骤**
1. 在 HPC 上准备数据
2. 创建公开配置: `configs/datasets/{dataset_id}.yaml`
3. 创建私有配置: `configs_private/datasets/{dataset_id}.yaml`
4. 创建 DE 对比: `configs/annotation/de_contrasts_{project}.yaml`（如需要）
5. 创建 Snakemake 入口: `workflow/datasets/{short_name}/`
6. 注册数据集: 添加到 `metadata/registry/dataset_registry.csv`
7. 运行流程: `cd workflow/datasets/{short_name}/ && snakemake --cores 16`

**详细说明**: 见 `docs/SOP_DATASET_ONBOARDING.md`

---

## **8. 调试流程失败的规则**

### **常见失败模式**

#### **1. 缺少输入文件**
```
Error: Input not found: results/scrna/novogene_llc/05_cluster/objects/seurat_clustered.rds
```
**解决方案**: 检查上一阶段是否完成，如需要则重新运行。

#### **2. 环境变量未设置**
```
Error: DS_PREFIX not set
```
**解决方案**:
```bash
export DS_CONFIG="configs_private/datasets/{dataset_id}.yaml"
export DS_PREFIX="{short_name}"
```

#### **3. 配置优先级混淆**
**解决方案**: 用 `load_qc_params()` 调试，查看最终生效的配置。

#### **4. CNV 输出契约违反**
```
Error: undefined columns: scevan_label
```
**解决方案**: 检查 CNV 模块输出格式，验证必需列是否存在。

#### **5. 细胞类型映射不同步**
```
⚠️  Unmapped clusters: 27, 28, 29, 30
```
**解决方案**: 用新 cluster 更新 `configs/annotation/celltype_mapping.csv`。

**更多详情**: 见完整文档第 8 节。

---

## **9. 命名约定**

### **数据集 ID**
```
{来源}_{供应商/登录号}_{疾病}_{模态}_{描述}_{日期}_{版本}

示例:
- INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1
- PUB_GEO_LUAD_scRNA_GSE253718_v1
```

### **数据集短名称（ds_prefix）**
```
{项目}_{物种/模型}

示例:
- novogene_llc
- gse253718
```

### **元数据列**
- **核心**: `sample_id`, `group`, `nFeature_RNA`, `nCount_RNA`, `percent.mt`
- **细胞周期**: `S.Score`, `G2M.Score`, `Phase`, `CC.Difference`
- **聚类**: `seurat_clusters`, `clusters_res0.3` 等
- **注释**: `celltype_L1`, `celltype_L2`, `annotation_confidence`, `is_artifact`
- **CNV**: `scevan_label`, `scevan_subclone`

---

## **10. 已知模糊点 / 待办事项**

### **需要确认**
1. Phase 09-11 优先级和目标细胞类型
2. Meta 分析策略（跨数据集整合方法）
3. 模板化时间线
4. CNV 模块版本控制机制

### **已知问题**
1. 备份文件（`.bak_20260409`）应删除
2. 文档空缺: CNV 契约、配置优先级指南、函数库文档
3. 注册表中的测试数据集（9 个 `VEND_Test_*` 条目）
4. 大对象脆弱性（13.5 GB RDS，无备份/校验和）

### **未来增强**
1. 契约验证函数
2. 所有编排器中的环境检查
3. 手动执行的快速运行脚本
4. 配置链和 CNV 契约的自动化测试

---

## **11. 新协作者接入：最短路径**

### **第 1 天：理解大局（30 分钟）**
1. 阅读 `README.md`
2. 阅读本文档（第 1-5 节）
3. 探索 `results/scrna/novogene_llc/`

### **第 2 天：设置环境（1-2 小时）**
1. 阅读 `docs/SETUP.md`
2. 创建数据符号链接: `ln -s /data1/morgen/data_center/Lung_Cancer_2026 data`
3. 恢复 R 环境: `Rscript -e 'renv::restore()'`
4. 激活 Snakemake: `conda activate snakemake_env`

### **第 3 天：运行测试阶段（2-3 小时）**
1. 手动重新运行 Phase 08 DE
2. 将输出与现有结果比较
3. 验证图表

### **第 4 天：理解配置系统（1-2 小时）**
1. 阅读第 4 节（数据集配置约定）
2. 用 `load_qc_params()` 调试配置优先级
3. 比较项目 vs 数据集配置

### **第 5 天：添加新数据集（练习）**
1. 阅读第 7 节（添加新数据集的规则）
2. 用 GSE253718 练习（已存在）
3. 在一个样本上运行 Phase 02 QC

---

## **12. 快速参考：文件路径**

### **关键对象**
```
# 唯一真值来源
results/scrna/{ds_prefix}/06_annotate/objects/seurat_annotated_final.rds

# DE 结果
results/scrna/{ds_prefix}/08_de/deseq2/all_de_results_combined.csv
results/scrna/{ds_prefix}/08_de/enrichment/all_enrichment_combined.csv
```

### **关键配置文件**
```
# 数据集配置
configs/datasets/{dataset_id}.yaml                           # 公开
configs_private/datasets/{dataset_id}.yaml                   # 私有

# 分析参数
configs/params/scrna_qc_params.yaml                          # QC
configs/params/scrna_annotation_params.yaml                  # 小鼠注释
configs/params/scrna_annotation_params_human.yaml            # 人类注释

# 人工注释
configs/annotation/celltype_mapping.csv

# DE 对比
configs/annotation/de_contrasts.yaml                         # INT_Novogene
configs/annotation/de_contrasts_gse253718.yaml               # GSE253718
```

### **工具函数**
```
# Layer 1: 项目级
scripts/utils/utils_io.R
scripts/utils/utils_plotting.R
scripts/utils/utils_registry.R

# Layer 2: scRNA 专用
workflow/scrna/functions/io_scrna.R
workflow/scrna/functions/qc_utils.R
workflow/scrna/functions/annotation_utils.R
```

---

## **13. 联系与协作**

### **主要分析师**
- **姓名**: morgen
- **角色**: 流程架构师，主要分析师
- **关注**: scRNA-seq 流程开发，INT_Novogene 分析

### **协作协议**
1. 修改代码前：讨论设计
2. 添加数据集前：遵循第 7 节
3. 更改配置前：理解优先级链（第 4 节）
4. 重大更改后：更新本文档和 README.md

---

## **14. 常见陷阱**

### **陷阱 1：运行脚本时未设置环境变量**
```bash
# ❌ 错误
Rscript workflow/scrna/08_de/run_de_pipeline.R

# ✅ 正确
export DS_CONFIG="configs_private/datasets/INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1.yaml"
export DS_PREFIX="novogene_llc"
Rscript workflow/scrna/08_de/run_de_pipeline.R
```

### **陷阱 2：配置覆盖未识别**
数据集配置覆盖项目配置。检查数据集 YAML 中的 `qc_overrides`。

### **陷阱 3：忘记质量门控**
Snakemake 在门控处停止。审查并批准: `touch .gate_approved`

### **陷阱 4：细胞类型映射不同步**
重新聚类后，用新 cluster 更新 `celltype_mapping.csv`。

### **陷阱 5：CNV 模块格式变化**
如果 CNV 输出列变化，更新 `07_cnv/01_cnv_to_annotation.R`。

---

## **15. 术语表**

| **术语** | **定义** |
|----------|---------------|
| **ds_prefix** | 数据集短名称（如 `novogene_llc`） |
| **dataset_id** | 完整标识符（如 `INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1`） |
| **Phase** | 主要分析阶段（01-11） |
| **哨兵文件** | `.phase{N}_done` 标记 |
| **质量门控** | 人工审查检查点（`.gate_approved`） |
| **Layer 1/2/3** | 函数层次（项目/scRNA/阶段特定） |
| **MAD** | 中位绝对偏差（自适应 QC） |
| **Pseudobulk** | 样本级聚合用于 DE |
| **celltype_L1/L2** | 细胞类型层次（广泛/特定） |
| **CNV** | 拷贝数变异 |
| **契约** | 接口规范（输入/输出格式） |

---

## **16. 紧急联系与资源**

### **出问题时**

| **问题** | **首要行动** | **升级** |
|-----------|-----------------|----------------|
| 流程失败 | 检查第 8 节 | 查看 `logs/snakemake/` 中的日志 |
| 配置不工作 | 用 `load_qc_params()` 调试 | 检查第 4 节 |
| 对象损坏 | 检查文件大小，尝试 `readRDS()` | 恢复备份或重新运行 |
| 磁盘空间不足 | 检查 `du -sh results/` | 清理或扩展存储 |
| R 包缺失 | `renv::restore()` | 检查 `renv.lock` |
| Snakemake DAG 错误 | `snakemake -np` | 检查 Snakefile 语法 |

### **关键文档**
- `README.md` — 项目概览
- `docs/SETUP.md` — 环境设置
- `docs/SOP_DATASET_ONBOARDING.md` — 数据集接入
- `workflow/scrna/README.md` — 流程架构
- 本文档 — 工作记忆

---

## **17. 版本历史**

| **版本** | **日期** | **变更** | **作者** |
|-------------|----------|-------------|-----------|
| 1.0.0 | 2026-04-20 | 初始工作记忆创建 | AI Assistant + morgen |

---

**项目协作说明书 v1.0.0 结束**

---

## **相关文档**

- **设置与环境**: `docs/SETUP.md`
- **数据集接入**: `docs/SOP_DATASET_ONBOARDING.md`
- **流程架构**: `workflow/scrna/README.md`
- **项目概览**: `README.md`
