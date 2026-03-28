# =============================================================================
# File   : scripts/00_setup_base_env.R
# Purpose: Initialize a reproducible base R environment for the atlas project
# Author : Morgen
# Note   : This script is for the ROOT project environment, not module-specific envs
# =============================================================================

# -----------------------------------------------------------------------------
# 0) Global options
# -----------------------------------------------------------------------------
# Why:
# - digits = 10：打印数值时保留更高精度，便于检查统计结果与中间对象
# - warn = 1：即时打印 warning（推荐），不要全局 suppress warnings
#             因为安装依赖和编译阶段的 warning 常常是后续报错的前兆
options(
  digits = 10,
  warn = 1
)

# -----------------------------------------------------------------------------
# 1) Mirror settings (CRAN / Bioconductor)
# -----------------------------------------------------------------------------
# Why:
# - 长期项目最怕“镜像不一致（mirror inconsistency）”导致依赖冲突
# - 建议 CRAN 与 Bioconductor 尽量使用同一体系下更稳定的镜像
# - 之前遇到过 dependency conflict，优先推荐 USTC
#
# If you insist on Westlake, you can switch BIOC_MIRROR manually.
CRAN_MIRROR <- "https://mirrors.ustc.edu.cn/CRAN/"
BIOC_MIRROR <- "https://mirrors.ustc.edu.cn/bioc/"
# 备选（如需测试）：BIOC_MIRROR <- "https://mirrors.westlake.edu.cn/bioconductor"

options(
  repos = c(CRAN = CRAN_MIRROR),
  BioC_mirror = BIOC_MIRROR
)

# Why:
# - BiocManager::repositories() 会把 BioCsoft / BioCann / BioCexp 等仓库补齐
# - 否则只写 options(repos=...) 很容易把 Bioconductor 标准仓库覆盖掉
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = CRAN_MIRROR)
}
options(repos = BiocManager::repositories())

# -----------------------------------------------------------------------------
# 2) Parallel settings
# -----------------------------------------------------------------------------
# Why:
# - HPC 资源充足不代表安装包一定要开满线程
# - 高并行安装会增加内存占用与编译波动
# - 12 核是比较稳妥的折中值（balanced setting）
#
# Rule of thumb:
# - 常规环境安装：8~12
# - 大型编译 / 重依赖场景：12~16
# - 不建议上来直接拉满 24+
options(Ncpus = 12)

# -----------------------------------------------------------------------------
# 3) Initialize renv (project-local library)
# -----------------------------------------------------------------------------
# Why:
# - renv::init(bare = TRUE) 只创建最小环境骨架，不会扫描整个项目目录
# - bare = TRUE 对大型科研仓库尤其重要，可避免首次初始化时“扫描过多对象/依赖”
# - init() 会自动 activate()，无需再手动 renv::activate()
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = CRAN_MIRROR)
}

if (!file.exists("renv.lock")) {
  renv::init(bare = TRUE)
} else {
  message("[INFO] renv.lock already exists; skipping renv::init().")
}

# -----------------------------------------------------------------------------
# 4) Enable pak + renv linkage
# -----------------------------------------------------------------------------
# Why:
# - pak 用于更快的依赖解析（dependency resolution）与安装
# - renv.config.pak.enabled = TRUE 可让 renv 与 pak 更好协同
options(renv.config.pak.enabled = TRUE)

# 打印当前库路径（library paths），便于确认是否已进入项目本地环境
message("[INFO] Current .libPaths():")
print(.libPaths())

# -----------------------------------------------------------------------------
# 5) Install pak if missing
# -----------------------------------------------------------------------------
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak", repos = CRAN_MIRROR)
}

# -----------------------------------------------------------------------------
# 6) Optional diagnostics
# -----------------------------------------------------------------------------
# Why:
# - 某些包是否可见、仓库是否工作正常，适合在搭环境初期做快速诊断
# - 这里默认关闭；需要时再打开，避免每次运行都打印一堆无关信息
RUN_DIAGNOSTIC_SEARCH <- FALSE

if (isTRUE(RUN_DIAGNOSTIC_SEARCH)) {
  message("[INFO] Running pak diagnostic search for 'qs' ...")
  tryCatch(
    pak::pkg_search("qs"),
    error = function(e) message("[WARN] pak::pkg_search('qs') failed: ", e$message)
  )
}

# -----------------------------------------------------------------------------
# 7) Optional pak extra helpers
# -----------------------------------------------------------------------------
# Why:
# - pak_install_extra() 不是主环境安装的硬依赖
# - 如果 pak 版本支持，可执行；若不存在则直接跳过
if ("pak_install_extra" %in% getNamespaceExports("pak")) {
  tryCatch(
    {
      pak::pak_install_extra()
      options(pak.no_extra_messages = TRUE)
      message("[INFO] pak extra tools installed.")
    },
    error = function(e) {
      message("[WARN] pak::pak_install_extra() failed or was skipped: ", e$message)
    }
  )
}

# -----------------------------------------------------------------------------
# 8) Package sets
# -----------------------------------------------------------------------------
# Why:
# - 将包分组（package grouping）有助于环境治理
# - 以后新增模块时，可明确知道“基础层”和“组学分析层”分别是什么

# ---- 8.1 Base packages: project IO / path / table handling / fast serialization
# Note:
# - qs2 是 qs 的继任者（successor），更适合长期项目
# - 暂不纳入 arrow：其编译链较重，且你之前已遇到 cmake 版本问题
base_pkgs <- c(
  "here",
  "fs",
  "data.table",
  "qs2"
)

# ---- 8.2 Omics packages: plotting / annotation / enrichment / bulk RNA / scRNA
omix_pkgs <- c(
  # Visualization
  "ggplot2",
  "patchwork",
  "cowplot",
  "ComplexHeatmap",
  "pheatmap",
  
  # Annotation / ID mapping
  "AnnotationDbi",
  "org.Hs.eg.db",
  
  # Enrichment / pathway analysis
  "clusterProfiler",
  "msigdbr",
  "fgsea",
  "enrichplot",
  
  # Bulk RNA-seq differential analysis
  "DESeq2",
  "edgeR",
  "limma",
  "sva",
  
  # Single-cell RNA-seq
  "Seurat",
  "harmony",
  "SingleCellExperiment",
  "scater",
  "scran",
  "scDblFinder",
  
  # Reporting
  "rmarkdown",
  "knitr"
)

# -----------------------------------------------------------------------------
# 9) Installation policy
# -----------------------------------------------------------------------------
# Why:
# - 对长期项目，不建议每次运行 setup 脚本都 upgrade = TRUE
# - 否则 lockfile 会不断漂移（version drift），影响可复现性
# - 首次安装可以设 TRUE；环境稳定后建议 FALSE
UPGRADE_PKGS <- FALSE

message("[INFO] Installing base packages ...")
pak::pkg_install(base_pkgs, upgrade = UPGRADE_PKGS, ask = FALSE)

message("[INFO] Installing omics packages ...")
pak::pkg_install(omix_pkgs, upgrade = UPGRADE_PKGS, ask = FALSE)

# -----------------------------------------------------------------------------
# 10) Module-specific packages
# -----------------------------------------------------------------------------
# Why:
# - CellChat 不应放在主环境里默认安装
# - 它属于模块化分析环境（module-specific environment）
# - 后续请在 modules/cellchat/scripts/00_setup_cellchat_env.R 中单独安装
#
# Example (DO NOT run here):
# pak::pkg_install("jinworks/CellChat@*release", ask = FALSE)

# -----------------------------------------------------------------------------
# 11) Snapshot environment
# -----------------------------------------------------------------------------
# Why:
# - renv::snapshot() 会把当前项目依赖写入 renv.lock
# - 这是“环境说明书（environment manifest）”，必须纳入版本控制
message("[INFO] Snapshotting environment ...")
renv::snapshot(prompt = FALSE)

message("[INFO] Checking renv status ...")
renv::status()

# -----------------------------------------------------------------------------
# 12) Quick sanity checks
# -----------------------------------------------------------------------------
# Why:
# - 安装成功不等于可正常加载（load）
# - 这里做最小关键包检测，快速发现问题
check_pkgs <- c("DESeq2", "limma", "Seurat", "clusterProfiler", "ComplexHeatmap")

pkg_load_status <- vapply(
  check_pkgs,
  function(pkg) {
    suppressPackageStartupMessages(require(pkg, character.only = TRUE))
  },
  FUN.VALUE = logical(1)
)

message("[INFO] Package load check:")
print(pkg_load_status)

if (!all(pkg_load_status)) {
  warning("Some key packages failed to load. Please inspect pkg_load_status carefully.")
} else {
  message("[INFO] Base environment setup completed successfully.")
}