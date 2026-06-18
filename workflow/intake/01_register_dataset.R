# workflow/intake/01_register_dataset.R
#
# 职责：
#   将 configs/datasets/<dataset_id>.yaml 中的配置解析并写入
#   metadata/registry/dataset_registry.csv，实现"配置即元数据"。
#
# 核心设计原则：
#   1. 配置驱动：所有数据集元数据从 YAML 读取，不手动编辑 CSV
#   2. 幂等性：同一 YAML 重复运行结果一致（Upsert 语义）
#   3. 门控前置：所有数据集进入 QC 流程前必须已登记
#   4. dataset_id 一致性：YAML 文件名（去后缀）必须与内部 dataset_id 字段完全一致
#
# 提供函数：
#   register_dataset_from_yaml(config_path)   ← 单文件注册（主函数）
#   register_all_configs(config_dir, ...)     ← 批量注册目录下所有 YAML
#
# 依赖：
#   scripts/utils_registry.R（提供 upsert_registry_entry / parse_nullable / %||%）
#
# 运行方式：
#   交互式（RStudio Web，推荐）：
#     source("workflow/intake/01_register_dataset.R")
#     register_dataset_from_yaml("configs/datasets/VEND_xxx.yaml")
#   命令行（仅在 RStudio Web 关联的 Rscript 下使用）：
#     Rscript workflow/intake/01_register_dataset.R configs/datasets/VEND_xxx.yaml
#
# 版本：1.1
# 作者：Morgen

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(tibble)
  library(dplyr)
})

# 加载注册表 CRUD 工具层
# 提供：upsert_registry_entry / parse_nullable / %||% / REGISTRY_PATH / SCHEMA_PATH
source(here("scripts", "utils", "utils_registry.R"))

# ==============================================================================
# 常量 (Constants)
# ==============================================================================

# 数据集 YAML 配置文件存放目录
CONFIG_ROOT <- here("configs", "datasets")

# ==============================================================================
# Section 1：YAML 字段校验辅助 (Validation Helpers)
# ==============================================================================

# 校验 dataset_id 与文件名是否一致
# 背景：dataset_id 是跨系统的核心键，文件名与字段值不一致会导致追溯混乱
# 参数：
#   cfg_id     YAML 内 dataset_id 字段的值
#   config_path  YAML 文件路径（用于提取文件名）
.check_id_matches_filename <- function(cfg_id, config_path) {
  filename_id <- tools::file_path_sans_ext(basename(config_path))
  if (!identical(cfg_id, filename_id)) {
    stop(sprintf(
      "[Error] dataset_id 与文件名不一致，注册终止。\n  文件名: %s\n  dataset_id: %s\n请确保两者完全一致（见 SOP 第 7 节）。",
      filename_id, cfg_id
    ))
  }
}

# 校验 files 节点结构
# 规则：
#   - files 必须存在且为命名列表（至少一个对象）
#   - 每个子节点必须包含至少一个格式键（qs / rds / csv 等）
#   - 格式键对应的值必须是非空字符串（文件名）
# 参数：
#   files_node  cfg$files（YAML 解析后的列表）
#   dataset_id  用于错误信息定位
.check_files_node <- function(files_node, dataset_id) {

  if (is.null(files_node) || !is.list(files_node) || length(files_node) == 0) {
    stop(sprintf(
      "[Error] 数据集 '%s' 的 `files` 节点为空或缺失。\n请参考 configs/datasets/_template.yaml 填写文件清单。",
      dataset_id
    ))
  }

  allowed_formats <- c("qs", "qs2", "rds", "rda", "csv", "tsv",
                        "h5ad", "loom", "mtx", "txt", "xlsx")

  for (obj_name in names(files_node)) {
    obj <- files_node[[obj_name]]

    # 每个对象必须是命名列表
    if (!is.list(obj) || is.null(names(obj))) {
      stop(sprintf(
        "[Error] 数据集 '%s'，files.%s 结构错误：期望命名列表（格式键 → 文件名），实际得到: %s",
        dataset_id, obj_name, class(obj)
      ))
    }

    # 每个格式键的值必须是非空字符串
    for (fmt in names(obj)) {
      fname <- obj[[fmt]]
      if (!is.character(fname) || nchar(trimws(fname)) == 0) {
        stop(sprintf(
          "[Error] 数据集 '%s'，files.%s.%s 文件名不能为空。",
          dataset_id, obj_name, fmt
        ))
      }
      # 警告未知格式键（不强制阻断，允许扩展）
      if (!fmt %in% allowed_formats) {
        warning(sprintf(
          "[Warning] 数据集 '%s'，files.%s.%s 使用了未在模板中定义的格式键 '%s'。",
          dataset_id, obj_name, fmt, fmt
        ))
      }
    }
  }

  message(sprintf(
    "[Info] files 节点校验通过：%d 个数据对象（%s）",
    length(files_node),
    paste(names(files_node), collapse = ", ")
  ))
}

# ==============================================================================
# Section 2：核心注册函数 (Core Registration)
# ==============================================================================

# 解析单个 YAML 配置文件并 Upsert 到注册表
#
# 幂等性保证：同一 dataset_id 重复运行 = 全量覆盖（不会产生重复行）
# 字段映射：YAML 嵌套结构 → 注册表平面结构
#           paths.raw_root      → hpc_raw_root
#           paths.manifest_root → hpc_manifest_root
#           paths.qc_root       → hpc_qc_root
#           paths.stage_root    → hpc_stage_root
#
# 参数：
#   config_path  YAML 配置文件路径（绝对路径或相对于项目根目录的路径）
# 返回：
#   写入注册表的 entry 列表（invisible）
register_dataset_from_yaml <- function(config_path,
                                      registry_path = REGISTRY_PATH) {

  # --- 步骤 1：路径存在性检查 ---
  config_path <- normalizePath(config_path, mustWork = FALSE)
  if (!file.exists(config_path)) {
    stop(sprintf("[Error] 配置文件不存在: %s", config_path))
  }

  message(sprintf(
    "\n[Process] ── 开始解析配置文件: %s ──",
    basename(config_path)
  ))

  # --- 步骤 2：YAML 解析 ---
  cfg <- tryCatch(
    yaml::read_yaml(config_path),
    error = function(e) {
      stop("[Error] YAML 格式错误（缩进/特殊字符问题）: ", conditionMessage(e))
    }
  )

  # --- 步骤 3：必填字段完整性校验 ---
  # 与 registry_schema.yaml 中 required = true 的字段对应
  required_top_level <- c(
    "dataset_id", "display_name", "dataset_type",
    "disease", "subtype", "modality", "data_level",
    "source_name", "owner",
    "status", "intake_qc_status", "harmonization_status", "analysis_qc_status",
    "register_date", "last_update",
    "paths", "files"
  )

  missing_fields <- setdiff(required_top_level, names(cfg))
  if (length(missing_fields) > 0) {
    stop(sprintf(
      "[Error] YAML 缺失以下必填字段: %s\n请参考 configs/datasets/_template.yaml 补充。",
      paste(missing_fields, collapse = ", ")
    ))
  }

  # --- 步骤 4：dataset_id 与文件名一致性校验 ---
  .check_id_matches_filename(cfg$dataset_id, config_path)

  # --- 步骤 5：paths 子节点存在性检查 ---
  if (is.null(cfg$paths$raw_root) || nchar(trimws(cfg$paths$raw_root)) == 0) {
    stop(sprintf(
      "[Error] 数据集 '%s' 的 paths.raw_root 为空，注册终止。",
      cfg$dataset_id
    ))
  }

  # --- 步骤 6：files 节点结构校验 ---
  .check_files_node(cfg$files, cfg$dataset_id)

  # --- 步骤 7：字段映射（YAML 嵌套 → 注册表平面结构）---
  #
  # 数值型字段使用 parse_nullable()，防止 YAML null/"NA" 写入字面量
  # 文本型可选字段使用 %||% 提供默认值
  entry <- list(

    # ── 身份标识 ──────────────────────────────────────────────────────────────
    dataset_id           = cfg$dataset_id,
    display_name         = cfg$display_name,
    dataset_type         = cfg$dataset_type,
    disease              = cfg$disease,
    subtype              = cfg$subtype,
    modality             = cfg$modality,
    data_level           = cfg$data_level,

    # ── 来源与归属 ────────────────────────────────────────────────────────────
    source_name          = cfg$source_name,
    source_accession     = parse_nullable(cfg$source_accession),
    owner                = cfg$owner,

    # ── 样本量声明（初始注册时允许为空）─────────────────────────────────────
    n_cohorts_claimed    = parse_nullable(cfg$n_cohorts_claimed),
    n_samples_claimed    = parse_nullable(cfg$n_samples_claimed),
    n_samples_aligned    = parse_nullable(cfg$n_samples_aligned),

    # ── 流程状态 ──────────────────────────────────────────────────────────────
    status               = cfg$status,
    intake_qc_status     = cfg$intake_qc_status,
    harmonization_status = cfg$harmonization_status,
    analysis_qc_status   = cfg$analysis_qc_status,

    # ── 分析角色与纳入决策 ────────────────────────────────────────────────────
    atlas_role           = cfg$atlas_role       %||% "undecided",
    inclusion_decision   = cfg$inclusion_decision %||% "pending",

    # ── 时间戳 ────────────────────────────────────────────────────────────────
    download_date        = parse_nullable(cfg$download_date),
    register_date        = cfg$register_date,
    last_update          = cfg$last_update,

    # ── HPC 路径（从 paths 嵌套节点展开为平面字段）──────────────────────────
    # 设计说明：注册表保持平面结构（便于 CSV 直接读写），
    #           YAML 保持嵌套结构（便于人类阅读和维护）
    hpc_raw_root         = cfg$paths$raw_root,
    hpc_manifest_root    = parse_nullable(cfg$paths$manifest_root),
    hpc_qc_root          = parse_nullable(cfg$paths$qc_root),
    hpc_stage_root       = parse_nullable(cfg$paths$stage_root),

    # ── 备注 ──────────────────────────────────────────────────────────────────
    notes                = parse_nullable(cfg$notes)
  )

 
  # ★ 唯一改动：把 registry_path 转发给 upsert
  upsert_registry_entry(entry, path = registry_path)

  message(sprintf("[OK] 数据集 '%s' 登记成功！", cfg$dataset_id))
  invisible(entry)
}

# ==============================================================================
# Section 3：批量注册函数 (Batch Registration)
# ==============================================================================

# 批量扫描并注册 configs/datasets/ 目录下所有 YAML 文件
#
# 排除规则：
#   - _ 开头的文件（如 _template.yaml）：模板/保留文件，不注册
#
# 参数：
#   config_dir    YAML 所在目录，默认 CONFIG_ROOT
#   pattern       文件名匹配正则，默认 "\\.yaml$"
#   stop_on_error 遇到单个文件报错时是否终止全局（默认 FALSE = 跳过并记录）
# 返回：
#   list，包含 success（成功 dataset_id 向量）和 failed（失败文件路径向量）
register_all_configs <- function(config_dir    = CONFIG_ROOT,
                                 pattern       = "\\.yaml$",
                                 stop_on_error = FALSE) {

  yaml_files <- list.files(config_dir, pattern = pattern, full.names = TRUE)

  # 排除 _ 开头的模板或保留文件（如 _template.yaml）
  yaml_files <- yaml_files[!grepl("/_[^/]+\\.yaml$", yaml_files)]

  if (length(yaml_files) == 0) {
    message(sprintf("[Warning] 目录 %s 下未找到可注册的 YAML 文件。", config_dir))
    return(invisible(list(success = character(0), failed = character(0))))
  }

  message(sprintf(
    "\n[Batch] 发现 %d 个 YAML 配置文件，开始批量注册...\n",
    length(yaml_files)
  ))

  success <- character(0)
  failed  <- character(0)

  for (f in yaml_files) {
    result <- tryCatch(
      {
        register_dataset_from_yaml(f)
        "ok"
      },
      error = function(e) {
        message(sprintf("[SKIP] %s 注册失败: %s", basename(f), conditionMessage(e)))
        "error"
      }
    )

    if (result == "ok") {
      success <- c(success, basename(f))
    } else {
      failed <- c(failed, f)
      if (stop_on_error) {
        stop(sprintf("[Fatal] 批量注册在 %s 处终止（stop_on_error = TRUE）。", basename(f)))
      }
    }
  }

  # 汇总报告
  message(sprintf(
    "\n[Batch Done] 成功: %d  失败: %d / 共: %d",
    length(success), length(failed), length(yaml_files)
  ))

  if (length(failed) > 0) {
    message("[Failed files]:\n  ", paste(basename(failed), collapse = "\n  "))
  }

  invisible(list(success = success, failed = failed))
}

# ==============================================================================
# Section 4：执行入口 (Execution Entry)
# ==============================================================================

if (interactive()) {
  # ── 交互式模式（RStudio Web，推荐）──────────────────────────────────────────
  # 取消注释对应行使用：
  #
  # 注册单个数据集：
  # register_dataset_from_yaml(
  #   here("configs", "datasets", "EXAMPLE_Lung_Bulk_PublicCohort_20260101_v1.yaml")
  # )
  #
  # 批量注册所有数据集：
  # register_all_configs()
  #
  # 弹窗选择单个文件（本地 Mac 调试用）：
  # register_dataset_from_yaml(file.choose())

  message(paste(
    "[Ready] 01_register_dataset.R 已加载。",
    "可用函数：",
    "  register_dataset_from_yaml('<path_to_yaml>')  ← 单文件注册",
    "  register_all_configs()                        ← 批量注册",
    sep = "\n"
  ))

} else if (sys.nframe() == 0) {
  # ── 命令行模式（仅在 RStudio Web 关联的 Rscript 下使用）────────────────────
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) == 0) {
    # 无参数：批量注册模式
    message("[CLI] 未指定文件路径，执行批量注册模式...")
    register_all_configs(stop_on_error = FALSE)

  } else if (length(args) == 1) {
    # 单参数：注册指定文件
    register_dataset_from_yaml(args[[1]])

  } else {
    stop(paste(
      "用法:",
      "  Rscript workflow/intake/01_register_dataset.R                      # 批量注册",
      "  Rscript workflow/intake/01_register_dataset.R <path_to_yaml>       # 单文件注册",
      sep = "\n"
    ))
  }
}
