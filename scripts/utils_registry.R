#' @title Dataset Registry Management Utilities
#' @description 提供注册表 CRUD 操作，确保食管癌/泛癌项目的元数据完整性。
#' @environment Ubuntu HPC / MacAir M4
#' @author Gemini Collaborative

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tibble)
  library(purrr)
})

# 配置项：注册表默认路径
REGISTRY_PATH <- here("metadata", "registry", "dataset_registry.csv")

# ==============================================================================
# 1. 基础 I/O 函数
# ==============================================================================

read_registry <- function(path = REGISTRY_PATH) {
  #' @description 读取数据集注册表 (Read Registry)
  #' @param path 注册表 CSV 路径
  #' @return A tibble.
  
  if (!file.exists(path)) {
    stop(sprintf("[Error] 注册表文件不存在: %s. 请检查路径或初始化文件。", path))
  }
  
  # 使用 col_types = cols(.default = "c") 强制按字符串读取，防止 ID 丢失前导零
  readr::read_csv(path, show_col_types = FALSE)
}

write_registry <- function(registry_df, path = REGISTRY_PATH) {
  #' @description 将注册表写回磁盘，并强制排序 (Write & Standardize)
  #' @param registry_df 待写入的 tibble/data.frame
  #' @param path 输出路径
  
  stopifnot(is.data.frame(registry_df))
  
  # 规范化：按 dataset_id 排序，确保 Git Diff 最小化
  registry_df <- registry_df %>% 
    dplyr::arrange(dataset_id)
  
  readr::write_csv(registry_df, path)
  
  # QC 清单 (QC Checklist):
  # 1. 重复检查: nrow(registry_df) == length(unique(registry_df$dataset_id))
  # 2. 缺失值检查: 确认关键字段没有空值
  # 3. 维度检查: 打印当前注册表规模
  message(sprintf("[Success] Registry updated: %d datasets recorded.", nrow(registry_df)))
  
  invisible(registry_df)
}

# ==============================================================================
# 2. 核心操作函数 (Upsert & Update)
# ==============================================================================

upsert_registry_entry <- function(entry, path = REGISTRY_PATH) {
  #' @description 插入或替换整行记录 (Insert or Replace Entry)
  #' @param entry 命名列表 (Named List) 或单行 tibble，必须包含 dataset_id
  #' @param path 注册表路径
  
  # 输入预处理 (Standardization)
  if (is.list(entry) && !inherits(entry, "data.frame")) {
    entry <- tibble::as_tibble(entry)
  }
  
  # 防御性校验 (Defensive Checks)
  if (!is.data.frame(entry) || nrow(entry) != 1) {
    stop("`entry` 必须是命名的 List 或单行 data.frame/tibble。")
  }
  
  if (!"dataset_id" %in% names(entry)) {
    stop("`entry` 缺失必需字段: `dataset_id`。")
  }
  
  registry_df <- read_registry(path)
  target_id <- as.character(entry$dataset_id[[1]])
  
  # 模式对齐 (Schema Alignment)
  # 1. 补充缺失列 (Fill missing columns with NA)
  missing_in_entry <- setdiff(names(registry_df), names(entry))
  for (col in missing_in_entry) {
    entry[[col]] <- NA
  }
  
  # 2. 丢弃多余列并发出警告 (Drop extra columns)
  extra_in_entry <- setdiff(names(entry), names(registry_df))
  if (length(extra_in_entry) > 0) {
    warning("丢弃了注册表中不存在的列: ", paste(extra_in_entry, collapse = ", "))
    entry <- entry %>% dplyr::select(all_of(names(registry_df)))
  }
  
  # 3. 确保列顺序完全一致
  entry <- entry[, names(registry_df), drop = FALSE]
  
  # 执行更新或插入 (Upsert Logic)
  # 移除旧 ID 的行 (如果存在)
  registry_df <- registry_df %>% 
    dplyr::filter(dataset_id != !!target_id)
  
  # 合并并写入
  registry_df <- dplyr::bind_rows(registry_df, entry)
  write_registry(registry_df, path)
}

update_registry_fields <- function(target_id, fields, path = REGISTRY_PATH) {
  #' @description 精确更新现有记录的特定字段 (Patch/Update Fields)
  #' @param target_id 目标 dataset_id (String)
  #' @param fields 包含更新内容的命名列表，例如 list(qc_status = "pass", samples_n = 48)
  #' @param path 注册表路径
  
  if (!is.list(fields) || length(fields) == 0) {
    stop("`fields` 必须是一个非空的命名列表。")
  }
  
  registry_df <- read_registry(path)
  
  # 定位索引 (Indexing)
  idx <- which(registry_df$dataset_id == target_id)
  if (length(idx) == 0) {
    stop(sprintf("[Error] 在注册表中未找到 ID: %s", target_id))
  }
  
  # 校验更新的列是否存在 (Schema Check)
  unknown_cols <- setdiff(names(fields), names(registry_df))
  if (length(unknown_cols) > 0) {
    stop("注册表中不存在以下列，更新失败: ", paste(unknown_cols, collapse = ", "))
  }
  
  # 批量更新字段
  for (nm in names(fields)) {
    registry_df[idx, nm] <- fields[[nm]]
  }
  
  write_registry(registry_df, path)
}

# ==============================================================================
# 3. 辅助函数
# ==============================================================================

view_registry <- function(path = REGISTRY_PATH) {
  #' @description 在 RStudio 中预览注册表
  registry_df <- read_registry(path)
  View(registry_df)
  invisible(registry_df)
}