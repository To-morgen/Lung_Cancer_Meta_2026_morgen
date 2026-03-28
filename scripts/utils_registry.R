# scripts/utils_registry.R
#
# 职责：数据集注册表（dataset_registry.csv）的 CRUD 工具层
#
# 提供：
#   - 读/写注册表（read_registry / write_registry）
#   - Schema 校验（validate_registry，对照 registry_schema.yaml）
#   - 记录插入与整行替换（upsert_registry_entry）
#   - 字段精确更新（update_registry_fields）
#   - 安全补丁（safe_patch：存在则 patch，不存在则 upsert）
#   - 查询辅助（get_registry_entry / view_registry）
#
# 被调用方：
#   workflow/intake/01_register_dataset.R
#   workflow/qc/01_vendor_curated_qc.R
#   （未来其他工作流脚本）
#
# 版本变更：
#   v1.0 → 初始版本
#   v1.1 → 引入 filelock 并发保护 / col_types 强制字符 /
#           validate_registry 对接 registry_schema.yaml /
#           新增 safe_patch / 修复 sprintf 参数错位 Bug /
#           修复 na_mask NA 传染问题 / 清理孤立代码块
#
# 运行环境：Ubuntu HPC (rlab) / MacAir M4
# 作者：Morgen
# 版本：1.1

suppressPackageStartupMessages({
  library(here)
  library(yaml)      # 用于 .load_schema() 读取 registry_schema.yaml
  library(readr)
  library(dplyr)
  library(tibble)
  library(purrr)
})

# FIX v1.1: 按需加载 filelock，不强制依赖
# 单机环境可接受无 filelock；HPC 多 Job 并发写入时强烈建议安装
HAS_FILELOCK <- requireNamespace("filelock", quietly = TRUE)

# ==============================================================================
# 常量定义 (Constants)
# ==============================================================================

# 注册表 CSV 文件路径（所有函数默认操作此文件）
REGISTRY_PATH <- here("metadata", "registry", "dataset_registry.csv")

# 注册表 Schema 定义文件路径
# 设计决策：schema 放入 configs/registry/ 而非 metadata/registry/
# 原因：configs/ 存放"规则"，metadata/ 存放"状态数据"，职责分离
SCHEMA_PATH <- here("configs", "registry", "registry_schema.yaml")

# ==============================================================================
# Section 0：内部工具函数 (Internal Helpers)
# ==============================================================================

# NULL / 空值合并运算符
# 用途：替代 rlang::%||%，消除对 rlang 的隐式依赖
# 语义：x 为 NULL、字符串 "NA"、或空字符串时，返回 y；否则返回 x
# FIX v1.1: 原代码使用 rlang::%||% 但未在 Imports 中声明
`%||%` <- function(x, y) {
  if (is.null(x) || identical(x, "NA") || identical(x, "")) y else x
}

# YAML NA 值标准化函数
# 背景：yaml::read_yaml() 将 YAML 裸值 NA 解析为字符串 "NA"，
#       而非 R 的 NA_character_，需要显式转换
# FIX v1.1: 原代码未处理此问题，导致注册表中出现字面量 "NA" 字符串
parse_nullable <- function(x, na_strings = c("NA", "na", "N/A", "")) {
  if (is.null(x) || x %in% na_strings) return(NA_character_)
  as.character(x)
}

# ==============================================================================
# Section 1：Schema 加载 (Schema Loader)
# ==============================================================================

# 内部函数：加载并解析 registry_schema.yaml
# 返回值：解析后的 schema 列表；schema 文件不存在或解析失败时返回 NULL
# 设计原则：失败时优雅降级，不阻断注册流程（仅跳过 schema 级校验）
.load_schema <- function(path = SCHEMA_PATH) {

  # FIX v1.1: 原版 sprintf 参数错位 Bug（path 作为第三个参数被忽略）
  if (!file.exists(path)) {
    warning(sprintf(
      "[Warning] Schema 文件不存在: %s\n跳过 schema 校验，仅执行基础检查。",
      path
    ))
    return(NULL)
  }

  tryCatch(
    yaml::read_yaml(path),
    error = function(e) {
      warning(
        "[Warning] Schema 文件解析失败，跳过 schema 校验: ",
        conditionMessage(e)
      )
      NULL
    }
  )
}

# ==============================================================================
# Section 2：注册表校验 (Registry Validation)
# ==============================================================================

# 注册表完整性校验函数
# 被 write_registry() 在写入前自动调用，也可手动调用做健康检查
#
# 校验层级（由轻到重）：
#   Level 0 — 维度信息（仅打印，不报错）
#   Level 1 — 重复 dataset_id（ERROR，硬性阻断写入）
#   Level 2 — Schema 列覆盖度（WARNING，schema 有但 CSV 没有的列）
#   Level 3 — required = true 字段缺失值（ERROR，硬性阻断）
#   Level 4 — required = warn 字段缺失值（WARNING，软性提示）
#   Level 5 — enum 字段非法值（ERROR）
#   Level 6 — date 字段格式校验（ERROR）
#   Level 7 — string pattern 校验，如 dataset_id 命名规范（ERROR）
#
# 参数：
#   registry_df  待校验的 tibble/data.frame
#   schema_path  schema YAML 文件路径
validate_registry <- function(registry_df, schema_path = SCHEMA_PATH) {

  stopifnot(is.data.frame(registry_df))

  # --- Level 0: 维度信息（仅供参考，不触发报错）---
  message(sprintf(
    "[Info] Registry dims: %d rows × %d cols",
    nrow(registry_df), ncol(registry_df)
  ))

  # --- Level 1: 重复 dataset_id（硬性阻断，任何时候不允许重复主键）---
  dup_ids <- registry_df$dataset_id[duplicated(registry_df$dataset_id)]
  if (length(dup_ids) > 0) {
    stop(sprintf(
      "[Error] 注册表存在重复 dataset_id，写入终止:\n  %s",
      paste(dup_ids, collapse = "\n  ")
    ))
  }

  # --- 加载 Schema 文件（失败则仅完成 Level 0-1 基础校验后返回）---
  schema <- .load_schema(schema_path)
  if (is.null(schema)) {
    message("[Info] 已完成基础校验（无 schema 文件，跳过字段级校验）。")
    return(invisible(TRUE))
  }

  col_defs      <- schema$columns
  schema_cols   <- names(col_defs)
  registry_cols <- names(registry_df)

  # --- Level 2: Schema 列覆盖度（schema 定义了但 CSV 中缺失的列）---
  # 常见原因：schema 版本升级后 CSV 尚未迁移
  missing_in_csv <- setdiff(schema_cols, registry_cols)
  if (length(missing_in_csv) > 0) {
    warning(sprintf(
      "[Warning] 以下 schema 列在注册表 CSV 中不存在（可能需要列迁移）:\n  %s",
      paste(missing_in_csv, collapse = "\n  ")
    ))
  }

  # 仅对 CSV 和 schema 均存在的列执行后续字段级校验
  cols_to_check <- intersect(schema_cols, registry_cols)

  errors   <- character(0)
  warnings <- character(0)

  for (col in cols_to_check) {
    def  <- col_defs[[col]]
    vals <- registry_df[[col]]
    req  <- def$required %||% "false"

    # FIX v1.1: 使用 %in% 替代 ==，避免 NA 传染导致 na_mask 出现 NA 值
    # 原版：is.na(vals) | vals == "" | vals == "NA"
    # 问题：vals == "" 在 vals 含 NA 时返回 NA，sum(na_mask) 结果不可靠
    na_mask   <- is.na(vals) | vals %in% c("", "NA")
    n_missing <- sum(na_mask)

    # --- Level 3 & 4: 必填性校验 ---
    if (n_missing > 0) {
      missing_ids <- registry_df$dataset_id[na_mask]
      msg <- sprintf(
        "字段 '%s' 存在 %d 个缺失值（dataset_id: %s）",
        col, n_missing, paste(missing_ids, collapse = ", ")
      )
      if (identical(req, TRUE) || identical(req, "true")) {
        errors <- c(errors, paste("[Error]", msg))       # Level 3
      } else if (identical(req, "warn")) {
        warnings <- c(warnings, paste("[Warning]", msg)) # Level 4
      }
    }

    # --- Level 5: enum 合法值校验 ---
    if (!is.null(def$type) && def$type == "enum" && !is.null(def$allowed_values)) {
      non_na_vals  <- vals[!na_mask]
      allowed      <- as.character(def$allowed_values)
      illegal_vals <- unique(non_na_vals[!non_na_vals %in% allowed])
      if (length(illegal_vals) > 0) {
        errors <- c(errors, sprintf(
          "[Error] 字段 '%s' 存在非法值: {%s}。允许值: {%s}",
          col,
          paste(illegal_vals, collapse = ", "),
          paste(allowed,      collapse = ", ")
        ))
      }
    }

    # --- Level 6: date 格式校验（期望格式 YYYY-MM-DD）---
    if (!is.null(def$type) && def$type == "date") {
      date_pattern <- def$pattern %||% "^\\d{4}-\\d{2}-\\d{2}$"
      non_na_vals  <- vals[!na_mask]
      bad_dates    <- non_na_vals[!grepl(date_pattern, non_na_vals)]
      if (length(bad_dates) > 0) {
        errors <- c(errors, sprintf(
          "[Error] 字段 '%s' 存在不符合 YYYY-MM-DD 格式的值: {%s}",
          col, paste(unique(bad_dates), collapse = ", ")
        ))
      }
    }

    # --- Level 7: string pattern 校验（如 dataset_id 命名规范）---
    # 仅当 schema 中定义了 pattern 且类型为 string 时触发
    if (!is.null(def$pattern) && (is.null(def$type) || def$type == "string")) {
      non_na_vals <- vals[!na_mask]
      bad_vals    <- non_na_vals[!grepl(def$pattern, non_na_vals, perl = TRUE)]
      if (length(bad_vals) > 0) {
        errors <- c(errors, sprintf(
          "[Error] 字段 '%s' 不符合命名规范 (pattern: %s):\n    {%s}",
          col, def$pattern, paste(unique(bad_vals), collapse = ", ")
        ))
      }
    }
  }

  # --- 汇总输出（先打 warnings，再决定是否 stop）---
  if (length(warnings) > 0) purrr::walk(warnings, message)

  if (length(errors) > 0) {
    stop(
      "[Validation Failed] 注册表校验发现以下错误，写入终止:\n",
      paste(errors, collapse = "\n")
    )
  }

  message(sprintf(
    "[OK] Schema 校验通过（检查了 %d 个字段，%d 个 warnings，0 个 errors）",
    length(cols_to_check), length(warnings)
  ))
  invisible(TRUE)
}

# ==============================================================================
# Section 3：基础 I/O 函数 (Read / Write)
# ==============================================================================

# 读取注册表 CSV
# 返回：全列强制字符串类型的 tibble（防止 readr 自动类型推断导致前导零丢失等问题）
# FIX v1.1: 原版仅 show_col_types = FALSE，类型推断仍在发生
read_registry <- function(path = REGISTRY_PATH) {

  if (!file.exists(path)) {
    stop(sprintf(
      "[Error] 注册表文件不存在: %s\n请检查路径或确认文件已初始化。",
      path
    ))
  }

  readr::read_csv(
    path,
    col_types      = readr::cols(.default = "c"),  # 强制所有列为字符串
    show_col_types = FALSE
  )
}

# 将注册表写回磁盘
# 写入前自动执行：validate_registry() → 排序（确保 Git diff 最小化）
# FIX v1.1: 原版写前校验仅在注释中提及，未实际调用
write_registry <- function(registry_df, path = REGISTRY_PATH) {

  stopifnot(is.data.frame(registry_df))

  # 写前校验（validate_registry 内部会对照 schema 进行完整性检查）
  validate_registry(registry_df)

  # 按 dataset_id 排序，保证每次写出结果确定性，减少 Git diff 噪声
  registry_df <- dplyr::arrange(registry_df, dataset_id)

  readr::write_csv(registry_df, path)
  message(sprintf(
    "[Success] Registry updated: %d datasets × %d fields → %s",
    nrow(registry_df), ncol(registry_df), basename(path)
  ))

  invisible(registry_df)
}

# ==============================================================================
# Section 4：核心操作函数 (Upsert / Patch / Safe Patch)
# ==============================================================================

# 插入或整行替换记录（Upsert）
#
# 适用场景：
#   - 首次注册一个新数据集（Insert）
#   - 全量覆盖某条已有记录（Replace）
#
# 注意：未传入的字段将被填充为 NA。
# 如需只更新部分字段，请使用 update_registry_fields() 或 safe_patch()。
#
# 参数：
#   entry  命名 List 或单行 tibble，必须包含 dataset_id
#   path   注册表路径
upsert_registry_entry <- function(entry, path = REGISTRY_PATH) {

  # 统一入参类型为单行 tibble
  if (is.list(entry) && !inherits(entry, "data.frame")) {
    entry <- tibble::as_tibble(entry)
  }

  if (!is.data.frame(entry) || nrow(entry) != 1) {
    stop("`entry` 必须是命名的 List 或单行 data.frame/tibble。")
  }

  if (!"dataset_id" %in% names(entry)) {
    stop("`entry` 缺失必需字段: `dataset_id`。")
  }

  # 文件锁（HPC 多 Job 并发写入保护）
  # FIX v1.1: 原版无并发保护，多 Job 同时写入会导致静默数据丢失
  lock_handle <- .acquire_lock(path)
  on.exit(.release_lock(lock_handle), add = TRUE)

  registry_df <- read_registry(path)
  target_id   <- as.character(entry$dataset_id[[1]])

  # Schema 对齐（三步）：
  # 1. 补充注册表有但 entry 缺失的列（填 NA）
  missing_in_entry <- setdiff(names(registry_df), names(entry))
  for (col in missing_in_entry) entry[[col]] <- NA_character_

  # 2. 警告并丢弃 entry 中不在注册表 schema 内的多余列
  extra_in_entry <- setdiff(names(entry), names(registry_df))
  if (length(extra_in_entry) > 0) {
    warning(
      "[Warning] 以下字段不在注册表 schema 中，已忽略: ",
      paste(extra_in_entry, collapse = ", ")
    )
    entry <- dplyr::select(entry, dplyr::all_of(names(registry_df)))
  }

  # 3. 强制列顺序与注册表一致
  entry <- entry[, names(registry_df), drop = FALSE]

  # Upsert：先删除旧记录（若存在），再追加新记录
  registry_df <- dplyr::filter(registry_df, dataset_id != !!target_id)
  registry_df <- dplyr::bind_rows(registry_df, entry)

  write_registry(registry_df, path)
}


# 字段精确更新（Patch）
#
# 适用场景：流程状态流转，例如 QC 完成后更新 intake_qc_status。
# 若 target_id 不存在，直接报错（不自动新建）。
# 如需"不存在则新建"语义，请使用 safe_patch()。
#
# 参数：
#   target_id  目标 dataset_id（字符串）
#   fields     更新内容的命名列表，例如 list(intake_qc_status = "pass")
#   path       注册表路径
update_registry_fields <- function(target_id, fields, path = REGISTRY_PATH) {

  if (!is.list(fields) || length(fields) == 0) {
    stop("`fields` 必须是一个非空的命名列表。")
  }

  lock_handle <- .acquire_lock(path)
  on.exit(.release_lock(lock_handle), add = TRUE)

  registry_df <- read_registry(path)

  # 定位目标行，不存在则报错
  idx <- which(registry_df$dataset_id == target_id)
  if (length(idx) == 0) {
    # FIX v1.1: 修复 sprintf 参数错位 Bug（原版 target_id 被忽略）
    stop(sprintf(
      "[Error] 注册表中未找到 dataset_id: '%s'。\n如需新增记录，请先运行 01_register_dataset.R 或改用 safe_patch()。",
      target_id
    ))
  }

  # 拒绝更新不存在于注册表 schema 中的列（防止静默引入错误列）
  unknown_cols <- setdiff(names(fields), names(registry_df))
  if (length(unknown_cols) > 0) {
    stop(
      "[Error] 以下字段不在注册表 schema 中，更新终止: ",
      paste(unknown_cols, collapse = ", ")
    )
  }

  # 批量更新目标行的指定字段
  for (nm in names(fields)) {
    registry_df[idx, nm] <- as.character(fields[[nm]])
  }

  write_registry(registry_df, path)
}


# 安全补丁（Safe Patch）
#
# 语义：记录存在 → 精确 Patch 指定字段；记录不存在 → 以最小字段自动新建
# 适用场景：QC / 标准化脚本需要回写状态字段，但不确定该 ID 是否已注册
# FIX v1.1: 原版 QC 脚本直接调用 update_registry_fields，若 ID 未注册则崩溃
#
# 参数：
#   target_id  目标 dataset_id（字符串）
#   fields     更新内容的命名列表
#   path       注册表路径
safe_patch <- function(target_id, fields, path = REGISTRY_PATH) {

  registry_df <- read_registry(path)

  if (target_id %in% registry_df$dataset_id) {
    # 记录已存在：精确更新指定字段，其余字段保持不变
    update_registry_fields(target_id, fields, path)

  } else {
    # 记录不存在：发出警告后以最小条目自动新建
    # FIX v1.1: 修复 sprintf 参数错位 Bug（原版 target_id 被忽略）
    warning(sprintf(
      "[Warning] dataset_id '%s' 不在注册表中，将以最小字段自动新建。\n建议先运行 01_register_dataset.R 完成正式注册。",
      target_id
    ))
    min_entry <- c(list(dataset_id = target_id), fields)
    upsert_registry_entry(min_entry, path)
  }
}

# ==============================================================================
# Section 5：查询辅助函数 (Query Helpers)
# ==============================================================================

# 在 RStudio Viewer 中预览注册表
view_registry <- function(path = REGISTRY_PATH) {
  registry_df <- read_registry(path)
  View(registry_df)
  invisible(registry_df)
}

# 提取单条记录
# 返回：单行 tibble；若 dataset_id 不存在，返回 NULL 并打印提示
get_registry_entry <- function(target_id, path = REGISTRY_PATH) {
  registry_df <- read_registry(path)
  result      <- dplyr::filter(registry_df, dataset_id == target_id)

  if (nrow(result) == 0) {
    message(sprintf("[Info] 注册表中不存在 dataset_id: %s", target_id))
    return(NULL)
  }
  result
}

# ==============================================================================
# Section 6：内部文件锁管理 (Internal Lock Helpers)
# ==============================================================================

# 获取文件锁（内部函数，仅由 upsert / update 调用）
# 若 filelock 未安装，降级为无锁模式并发出警告
.acquire_lock <- function(path) {
  if (!HAS_FILELOCK) {
    warning(
      "[Warning] `filelock` 包未安装，跳过并发保护。",
      "单机环境可接受；HPC 环境建议: install.packages('filelock')"
    )
    return(NULL)
  }
  lock_path <- paste0(path, ".lock")
  filelock::lock(lock_path, timeout = 15000)  # 最长等待 15 秒
}

# 释放文件锁（内部函数，由 on.exit() 自动调用，确保异常时也能释放）
.release_lock <- function(lock_handle) {
  if (!is.null(lock_handle) && HAS_FILELOCK) {
    filelock::unlock(lock_handle)
  }
}
