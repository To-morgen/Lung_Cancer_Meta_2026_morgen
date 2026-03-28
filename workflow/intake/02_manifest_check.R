# workflow/intake/02_manifest_check.R
#
# 职责：
#   对照 configs/datasets/<dataset_id>.yaml 的 files 节点，
#   核查 HPC 上原始文件是否真实存在，并生成结构化文件清单（manifest CSV）。
#
# 核心设计原则：
#   1. 门控前置：canonical 格式文件缺失 → ERROR，阻断后续 QC 流程
#   2. 完整记录：manifest CSV 作为数据集文件状态的唯一真实来源
#   3. 幂等性：重复运行覆盖旧 manifest，不会产生重复记录
#   4. 注册表回写：检查完成后自动更新 hpc_manifest_root 字段
#
# 提供函数：
#   check_dataset_manifest(dataset_id_or_path)  ← 单数据集检查（主函数）
#   check_all_manifests(config_dir, ...)        ← 批量检查
#
# 依赖：
#   scripts/utils_registry.R
#   configs/datasets/<dataset_id>.yaml（提供 paths 和 files 声明）
#
# 输出：
#   <hpc_qc_root>/manifest/<dataset_id>_manifest.csv
#
# 版本：1.0
# 作者：Morgen

suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(readr)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(fs)       # 跨平台文件系统操作（file_size / file_info）
})

source(here("scripts", "utils_registry.R"))

# ==============================================================================
# 常量 (Constants)
# ==============================================================================

CONFIG_ROOT   <- here("configs", "datasets")

# manifest 默认输出子目录名（相对于 hpc_qc_root）
MANIFEST_SUBDIR <- "manifest"

# ==============================================================================
# Section 1：内部辅助函数 (Internal Helpers)
# ==============================================================================

# 根据 dataset_id 定位对应的 YAML 配置文件路径
# 若找不到则 stop
.resolve_config_path <- function(dataset_id, config_dir = CONFIG_ROOT) {
  path <- file.path(config_dir, paste0(dataset_id, ".yaml"))
  if (!file.exists(path)) {
    stop(sprintf(
      "[Error] 找不到数据集 '%s' 的配置文件: %s\n请先运行 01_register_dataset.R 完成注册。",
      dataset_id, path
    ))
  }
  path
}

# 从 YAML files 节点构建待检查的文件清单（tibble 格式，每行一个文件）
# 参数：
#   files_node  cfg$files（命名列表，object_name → list(format → filename)）
#   raw_root    数据集原始文件根目录（来自 cfg$paths$raw_root）
#   dataset_id  用于填充 dataset_id 列
# 返回：
#   tibble，含 object_name / format / is_canonical / filename / full_path
.build_file_table <- function(files_node, raw_root, dataset_id) {

  rows <- purrr::imap_dfr(files_node, function(obj, obj_name) {
    formats <- names(obj)
    purrr::imap_dfr(obj, function(fname, fmt) {
      tibble::tibble(
        dataset_id   = dataset_id,
        object_name  = obj_name,
        format       = fmt,
        # canonical = 每个 object 的第一个格式键（首选读取格式）
        is_canonical = identical(fmt, formats[[1]]),
        filename     = fname,
        full_path    = normalizePath(
          file.path(raw_root, fname),
          mustWork = FALSE
        )
      )
    })
  })

  rows
}

# 对 file_table 的每一行做文件存在性及属性检查
# 返回：
#   原 tibble 追加 exists / file_size_mb / last_modified / check_timestamp
.check_file_existence <- function(file_table) {

  check_ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  file_table %>%
    dplyr::mutate(
      exists = purrr::map_lgl(full_path, file.exists),

      # 文件大小（MB），不存在时为 NA
      file_size_mb = purrr::map_dbl(full_path, function(p) {
        if (!file.exists(p)) return(NA_real_)
        round(as.numeric(fs::file_size(p)) / 1024^2, 3)
      }),

      # 最后修改时间，不存在时为 NA
      last_modified = purrr::map_chr(full_path, function(p) {
        if (!file.exists(p)) return(NA_character_)
        format(fs::file_info(p)$modification_time, "%Y-%m-%d %H:%M:%S")
      }),

      check_timestamp = check_ts
    )
}

# ==============================================================================
# Section 2：Manifest 报告生成 (Report Generation)
# ==============================================================================

# 打印 manifest 检查摘要到控制台
# 包含：总文件数 / 存在数 / 缺失文件列表（按 canonical 优先展示）
.print_manifest_summary <- function(manifest_df, dataset_id) {

  n_total    <- nrow(manifest_df)
  n_exist    <- sum(manifest_df$exists)
  n_missing  <- n_total - n_exist

  message(sprintf(
    "\n[Manifest] %s — 文件核查结果: %d 存在 / %d 缺失 / %d 总计",
    dataset_id, n_exist, n_missing, n_total
  ))

  if (n_missing == 0) {
    message("[OK] 所有声明文件均已确认存在于 HPC。")
    return(invisible(NULL))
  }

  # 分 canonical 和非 canonical 分别列出
  missing_df       <- dplyr::filter(manifest_df, !exists)
  missing_canon    <- dplyr::filter(missing_df,  is_canonical)
  missing_noncanon <- dplyr::filter(missing_df, !is_canonical)

  if (nrow(missing_canon) > 0) {
    message("[ERROR] 以下 canonical（首选）文件缺失，将阻断 QC 流程:")
    purrr::walk(missing_canon$full_path, ~ message("    ✗ ", .x))
  }

  if (nrow(missing_noncanon) > 0) {
    message("[WARNING] 以下非 canonical（备用）文件缺失（不影响 QC）:")
    purrr::walk(missing_noncanon$full_path, ~ message("    △ ", .x))
  }

  invisible(NULL)
}

# ==============================================================================
# Section 3：核心检查函数 (Core Check)
# ==============================================================================

# 对单个数据集执行完整 manifest 检查
#
# 参数：
#   dataset_id_or_path  dataset_id 字符串，或 YAML 配置文件路径（两者均支持）
#   config_dir          YAML 配置目录，默认 CONFIG_ROOT
#   registry_path       注册表路径，默认 REGISTRY_PATH（来自 utils_registry.R）
#
# 返回：
#   manifest tibble（invisible）；canonical 文件缺失时 stop
#
# 副作用：
#   1. 在 <hpc_qc_root>/manifest/ 写出 <dataset_id>_manifest.csv
#   2. 回写注册表 hpc_manifest_root 字段
check_dataset_manifest <- function(dataset_id_or_path,
                                   config_dir    = CONFIG_ROOT,
                                   registry_path = REGISTRY_PATH) {

  # --- 步骤 1：解析入参，允许直接传 dataset_id 或 yaml 路径 ---
  if (file.exists(dataset_id_or_path)) {
    # 传入的是文件路径
    config_path <- normalizePath(dataset_id_or_path, mustWork = TRUE)
    dataset_id  <- tools::file_path_sans_ext(basename(config_path))
  } else {
    # 传入的是 dataset_id 字符串
    dataset_id  <- dataset_id_or_path
    config_path <- .resolve_config_path(dataset_id, config_dir)
  }

  message(sprintf(
    "\n[Manifest] ── 开始检查: %s ──", dataset_id
  ))

  # --- 步骤 2：读取 YAML ---
  cfg <- tryCatch(
    yaml::read_yaml(config_path),
    error = function(e) {
      stop("[Error] YAML 解析失败: ", conditionMessage(e))
    }
  )

  # --- 步骤 3：paths 与 files 节点基础校验 ---
  raw_root <- cfg$paths$raw_root
  if (is.null(raw_root) || nchar(trimws(raw_root)) == 0) {
    stop(sprintf("[Error] 数据集 '%s' 的 paths.raw_root 为空。", dataset_id))
  }

  # raw_root 目录本身是否存在
  if (!dir.exists(raw_root)) {
    stop(sprintf(
      "[Error] raw_root 目录不存在: %s\n请确认数据已下载到 HPC，或检查路径配置。",
      raw_root
    ))
  }

  if (is.null(cfg$files) || length(cfg$files) == 0) {
    stop(sprintf("[Error] 数据集 '%s' 的 files 节点为空，无法执行 manifest 检查。",
                 dataset_id))
  }

  # --- 步骤 4：构建文件清单表 ---
  file_table <- .build_file_table(cfg$files, raw_root, dataset_id)

  # --- 步骤 5：逐文件存在性检查 ---
  manifest_df <- .check_file_existence(file_table)

  # --- 步骤 6：打印摘要 ---
  .print_manifest_summary(manifest_df, dataset_id)

  # --- 步骤 7：canonical 文件缺失 → 硬性阻断 ---
  missing_canonical <- manifest_df %>%
    dplyr::filter(is_canonical, !exists)

  if (nrow(missing_canonical) > 0) {
    stop(sprintf(
      "[Error] 数据集 '%s' 存在 %d 个 canonical 文件缺失，manifest 检查终止。\n",
      "请先确认文件已完整下载，再重新运行 02_manifest_check.R。",
      dataset_id, nrow(missing_canonical)
    ))
  }

  # --- 步骤 8：写出 manifest CSV ---
  # 优先使用 YAML 中声明的 qc_root；若未声明则使用 raw_root 旁的 qc 子目录
  qc_root      <- cfg$paths$qc_root %||%
                  file.path(dirname(raw_root), "qc", dataset_id)
  manifest_dir <- file.path(qc_root, MANIFEST_SUBDIR)

  # 创建目录（若不存在）
  if (!dir.exists(manifest_dir)) {
    dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)
    message(sprintf("[Info] 创建 manifest 输出目录: %s", manifest_dir))
  }

  manifest_path <- file.path(
    manifest_dir,
    sprintf("%s_manifest.csv", dataset_id)
  )

  readr::write_csv(manifest_df, manifest_path)
  message(sprintf("[Info] Manifest 已写出: %s", manifest_path))

  # --- 步骤 9：回写 hpc_manifest_root 到注册表 ---
  # 使用 safe_patch：若注册表中已存在该 ID 则 patch，否则发出警告（不应发生）
  safe_patch(
    target_id = dataset_id,
    fields    = list(
      hpc_manifest_root = manifest_dir,
      last_update       = format(Sys.Date(), "%Y-%m-%d")
    ),
    path = registry_path
  )

  message(sprintf(
    "[OK] %s manifest 检查完成（%d 个文件，%d canonical，全部存在）。",
    dataset_id,
    nrow(manifest_df),
    sum(manifest_df$is_canonical)
  ))

  invisible(manifest_df)
}

# ==============================================================================
# Section 4：批量检查函数 (Batch Check)
# ==============================================================================

# 批量对 configs/datasets/ 目录下所有已注册数据集执行 manifest 检查
#
# 参数：
#   config_dir    YAML 所在目录，默认 CONFIG_ROOT
#   stop_on_error 单个数据集失败是否终止全局（默认 FALSE = 跳过并记录）
# 返回：
#   list(success = ..., failed = ...)
check_all_manifests <- function(config_dir    = CONFIG_ROOT,
                                stop_on_error = FALSE) {

  yaml_files <- list.files(config_dir, pattern = "\\.yaml$", full.names = TRUE)

  # 排除模板文件（_ 开头）
  yaml_files <- yaml_files[!grepl("/_[^/]+\\.yaml$", yaml_files)]

  if (length(yaml_files) == 0) {
    message("[Warning] 未找到可检查的 YAML 配置文件。")
    return(invisible(list(success = character(0), failed = character(0))))
  }

  message(sprintf(
    "\n[Batch Manifest] 发现 %d 个数据集，开始批量检查...\n",
    length(yaml_files)
  ))

  success <- character(0)
  failed  <- character(0)

  for (f in yaml_files) {
    dataset_id <- tools::file_path_sans_ext(basename(f))
    result <- tryCatch(
      {
        check_dataset_manifest(f)
        "ok"
      },
      error = function(e) {
        message(sprintf("[SKIP] %s 检查失败: %s", dataset_id, conditionMessage(e)))
        "error"
      }
    )

    if (result == "ok") {
      success <- c(success, dataset_id)
    } else {
      failed <- c(failed, dataset_id)
      if (stop_on_error) {
        stop(sprintf("[Fatal] 批量检查在 %s 处终止（stop_on_error = TRUE）。", dataset_id))
      }
    }
  }

  message(sprintf(
    "\n[Batch Done] 成功: %d  失败: %d / 共: %d",
    length(success), length(failed), length(yaml_files)
  ))

  if (length(failed) > 0) {
    message("[Failed]:\n  ", paste(failed, collapse = "\n  "))
  }

  invisible(list(success = success, failed = failed))
}

# ==============================================================================
# Section 5：执行入口 (Execution Entry)
# ==============================================================================

if (interactive()) {
  # ── 交互式模式（RStudio Web，推荐）──────────────────────────────────────────
  # 检查单个数据集（传 dataset_id 或 yaml 路径均可）：
  # check_dataset_manifest("VEND_AIshixin_LUAD_Bulk_22Cohorts_3019_20260307_v1")
  # check_dataset_manifest(here("configs", "datasets", "VEND_xxx.yaml"))
  #
  # 批量检查所有数据集：
  # check_all_manifests()

  message(paste(
    "[Ready] 02_manifest_check.R 已加载。",
    "可用函数：",
    "  check_dataset_manifest('<dataset_id 或 yaml 路径>')  ← 单数据集检查",
    "  check_all_manifests()                               ← 批量检查",
    sep = "\n"
  ))

} else if (sys.nframe() == 0) {
  # ── 命令行模式 ────────────────────────────────────────────────────────────
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) == 0) {
    # 无参数：批量检查模式
    message("[CLI] 未指定数据集，执行批量 manifest 检查...")
    check_all_manifests(stop_on_error = FALSE)

  } else if (length(args) == 1) {
    # 单参数：检查指定数据集
    check_dataset_manifest(args[[1]])

  } else {
    stop(paste(
      "用法:",
      "  Rscript workflow/intake/02_manifest_check.R                       # 批量检查",
      "  Rscript workflow/intake/02_manifest_check.R <dataset_id 或 path>  # 单数据集",
      sep = "\n"
    ))
  }
}
