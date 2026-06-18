# workflow/bulk/qc/01_vendor_curated_qc.R
#
# 职责：
#   对 vendor_curated 类数据集执行 Intake QC 门控检查（Gate 1）。
#   加载 canonical 格式的各数据对象，执行针对性质量检查，
#   生成 QC 报告 CSV，并回写 intake_qc_status 到注册表。
#
# 适用数据集类型：
#   dataset_type = vendor_curated
#   modality     = bulk
#
# 核心 Gate 规则：
#   fail    → 注册表写入 intake_qc_status = fail，阻断后续 harmonize 流程
#   warning → 注册表写入 intake_qc_status = warning，允许继续但需人工确认
#   pass    → 注册表写入 intake_qc_status = pass，可进入 harmonize
#
# 提供函数：
#   run_vendor_qc(dataset_id_or_path)     ← 单数据集 QC（主函数）
#   run_all_vendor_qc(config_dir, ...)    ← 批量 QC
#
# 依赖：
#   scripts/utils_registry.R
#   qs（加载 .qs 格式文件，须提前安装）
#
# 输出：
#   <hpc_qc_root>/intake/<dataset_id>_intake_qc.csv
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
  library(tidyr)
})

source(here("scripts", "utils", "utils_registry.R"))

# qs 包按需加载（加载 .qs 格式必须）
HAS_QS <- requireNamespace("qs", quietly = TRUE)

# ==============================================================================
# 常量 (Constants)
# ==============================================================================

CONFIG_ROOT <- here("configs", "datasets")

# QC 报告输出子目录（相对于 hpc_qc_root）
QC_SUBDIR <- "intake"

# 临床元数据中预期存在的关键字段（缺失时触发 warning）
# 按项目需求在此处集中维护，不要散落在函数内部
CLINICAL_KEY_FIELDS <- c(
  "sample_id",    # 样本唯一标识（必须）
  "cohort",       # 队列标识（必须）
  "OS",           # 总生存状态（0/1）
  "OS.time",      # 总生存时间（月/天）
  "gender",       # 性别
  "age",          # 年龄
  "stage"         # 肿瘤分期
)

# 表达矩阵中可能存放基因名的列名（按优先级排序，第一个命中者生效）
GENE_ID_CANDIDATES <- c("gene", "gene_id", "gene_name", "symbol", "Gene", "GeneSymbol")

# 样本数声明偏差容忍阈值（超出则触发 warning）
SAMPLE_COUNT_TOLERANCE <- 0.05   # 5%

# 临床关键字段缺失率容忍阈值（超出则触发 warning）
MISSING_RATE_THRESHOLD <- 0.20   # 20%

# ==============================================================================
# Section 1：内部工具函数 (Internal Helpers)
# ==============================================================================

# 初始化 QC 记录收集器（返回一个闭包对，用于追加记录和导出结果）
# 设计意图：统一收集 info / warning / error 三种级别的检查结果，
#           避免函数间反复传递 errors / warnings 向量
.make_qc_collector <- function(dataset_id) {

  records <- list()  # 内部状态：QC 记录列表

  # 追加一条 QC 记录
  add <- function(check_name, object_name, level,
                  detail, value = NA_character_) {
    records[[length(records) + 1]] <<- tibble::tibble(
      dataset_id      = dataset_id,
      check_name      = check_name,
      object_name     = object_name,
      level           = level,   # "info" / "warning" / "error"
      detail          = detail,
      value           = as.character(value),
      check_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
  }

  # 导出所有记录为 tibble
  export <- function() {
    if (length(records) == 0) return(tibble::tibble())
    dplyr::bind_rows(records)
  }

  # 判断当前最高严重等级
  get_status <- function() {
    df <- export()
    if (nrow(df) == 0 || !"error" %in% df$level) {
      if ("warning" %in% df$level) return("warning")
      return("pass")
    }
    "fail"
  }

  list(add = add, export = export, get_status = get_status)
}


# 根据 dataset_id 定位 YAML 路径（复用 manifest_check 的同名内部函数逻辑）
.resolve_config_path <- function(dataset_id, config_dir = CONFIG_ROOT) {
  path <- file.path(config_dir, paste0(dataset_id, ".yaml"))
  if (!file.exists(path)) {
    stop(sprintf(
      "[Error] 找不到数据集 '%s' 的配置文件: %s",
      dataset_id, path
    ))
  }
  path
}


# 加载单个数据对象文件（支持 qs / rds / csv / tsv）
# 参数：
#   full_path  文件绝对路径
#   format     格式键（"qs" / "rds" / "csv" / "tsv"）
# 返回：
#   加载后的 R 对象（data.frame / matrix / list 等）
.load_data_object <- function(full_path, format) {

  if (!file.exists(full_path)) {
    stop(sprintf("[Error] 文件不存在（manifest 检查应已拦截）: %s", full_path))
  }

  switch(format,
    "qs" = {
      if (!HAS_QS) stop("[Error] 加载 .qs 格式需要安装 `qs` 包: install.packages('qs')")
      qs::qread(full_path)
    },
    "qs2" = {
      if (!HAS_QS) stop("[Error] 加载 .qs2 格式需要安装 `qs` 包: install.packages('qs')")
      qs::qread(full_path)
    },
    "rds" = {
      readRDS(full_path)
    },
    "rda" = {
      e <- new.env(parent = emptyenv())
      load(full_path, envir = e)
      # .rda 可能含多个对象，取第一个并发出提示
      obj_names <- ls(e)
      if (length(obj_names) > 1) {
        warning(sprintf("[Warning] .rda 文件含 %d 个对象，仅加载第一个: %s",
                        length(obj_names), obj_names[[1]]))
      }
      get(obj_names[[1]], envir = e)
    },
    "csv" = {
      readr::read_csv(full_path, col_types = readr::cols(), show_col_types = FALSE)
    },
    "tsv" = {
      readr::read_tsv(full_path, col_types = readr::cols(), show_col_types = FALSE)
    },
    stop(sprintf("[Error] 不支持的格式: '%s'", format))
  )
}


# 识别表达矩阵中的基因 ID 列
# 参数：
#   df  表达矩阵（data.frame 格式，基因为行，但可能混有基因名列）
# 返回：
#   基因 ID 列名（字符串），若未找到则返回 NULL
.detect_gene_col <- function(df) {
  hit <- intersect(GENE_ID_CANDIDATES, names(df))
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

# ==============================================================================
# Section 2：各对象类型 QC 检查函数
# ==============================================================================

# ── 2a. 表达矩阵 QC（symbol / counts / tpm 等）─────────────────────────────

# 检查表达矩阵质量
# 参数：
#   expr_obj     加载后的表达矩阵（data.frame 或 matrix）
#   object_name  逻辑对象名（用于报告）
#   qc           QC 收集器（.make_qc_collector 返回的列表）
# 返回：
#   提取的样本 ID 向量（用于后续交叉比对）
.check_expression <- function(expr_obj, object_name, qc) {

  # ── 类型检查：统一转为 data.frame 操作 ──
  if (inherits(expr_obj, "matrix")) {
    expr_df <- as.data.frame(expr_obj)
  } else if (is.data.frame(expr_obj)) {
    expr_df <- expr_obj
  } else {
    qc$add(
      check_name  = "expression_type",
      object_name = object_name,
      level       = "error",
      detail      = sprintf("期望 data.frame 或 matrix，实际类型: %s", class(expr_obj)[[1]])
    )
    return(character(0))
  }

  # ── 维度信息 ──
  gene_col    <- .detect_gene_col(expr_df)
  numeric_cols <- names(expr_df)[sapply(expr_df, is.numeric)]
  n_genes     <- nrow(expr_df)
  n_samples   <- length(numeric_cols)

  qc$add(
    check_name  = "expression_dims",
    object_name = object_name,
    level       = "info",
    detail      = sprintf("基因数: %d，数值列（样本）数: %d", n_genes, n_samples),
    value       = sprintf("%d x %d", n_genes, n_samples)
  )

  # ── 基因 ID 列检测 ──
  if (is.null(gene_col)) {
    qc$add(
      check_name  = "expression_gene_col",
      object_name = object_name,
      level       = "warning",
      detail      = sprintf(
        "未找到标准基因 ID 列，候选列名: %s。请检查列命名规范。",
        paste(GENE_ID_CANDIDATES, collapse = " / ")
      )
    )
  } else {
    # 基因 ID 唯一性
    n_dup_genes <- sum(duplicated(expr_df[[gene_col]]))
    if (n_dup_genes > 0) {
      qc$add(
        check_name  = "expression_gene_uniqueness",
        object_name = object_name,
        level       = "error",
        detail      = sprintf("基因 ID 列 '%s' 存在 %d 个重复值", gene_col, n_dup_genes),
        value       = n_dup_genes
      )
    } else {
      qc$add(
        check_name  = "expression_gene_uniqueness",
        object_name = object_name,
        level       = "info",
        detail      = "基因 ID 唯一性：通过"
      )
    }
  }

  # ── NA / Inf 检查（仅数值列）──
  na_count  <- sum(is.na(expr_df[, numeric_cols, drop = FALSE]))
  inf_count <- sum(is.infinite(as.matrix(expr_df[, numeric_cols, drop = FALSE])))

  if (na_count > 0) {
    qc$add(
      check_name  = "expression_na",
      object_name = object_name,
      level       = "error",
      detail      = sprintf("表达矩阵含 %d 个 NA 值", na_count),
      value       = na_count
    )
  } else {
    qc$add(
      check_name  = "expression_na",
      object_name = object_name,
      level       = "info",
      detail      = "NA 检查：无 NA 值"
    )
  }

  if (inf_count > 0) {
    qc$add(
      check_name  = "expression_inf",
      object_name = object_name,
      level       = "error",
      detail      = sprintf("表达矩阵含 %d 个 Inf/-Inf 值", inf_count),
      value       = inf_count
    )
  }

  # ── 负值检查 ──
  n_neg <- sum(expr_df[, numeric_cols, drop = FALSE] < 0, na.rm = TRUE)
  if (n_neg > 0) {
    qc$add(
      check_name  = "expression_negative",
      object_name = object_name,
      level       = "warning",
      detail      = sprintf(
        "检测到 %d 个负值。若为 TPM/counts 数据，负值不符合预期；",
        "若为 log2 fold change 数据，可忽略。",
        n_neg
      ),
      value       = n_neg
    )
  }

  # ── 尺度检测：判断是否已 log 转换 ──
  # 启发式规则：若所有数值列的中位数 < 30，大概率已 log 转换
  col_medians <- sapply(expr_df[, numeric_cols, drop = FALSE],
                        function(x) median(x, na.rm = TRUE))
  overall_median <- median(col_medians, na.rm = TRUE)
  scale_hint <- if (overall_median < 30) "可能已 log 转换" else "可能为原始计数/TPM（未 log）"
  qc$add(
    check_name  = "expression_scale_hint",
    object_name = object_name,
    level       = "info",
    detail      = sprintf("中位数 = %.2f，%s（请结合供应商说明确认）", overall_median, scale_hint),
    value       = round(overall_median, 2)
  )

  # 返回样本 ID（数值列名）供后续交叉比对
  numeric_cols
}


# ── 2b. 临床元数据 QC ───────────────────────────────────────────────────────

# 检查临床元数据质量
# 参数：
#   clinical_obj   加载后的临床数据（data.frame / tibble）
#   cfg            完整 YAML 配置（用于获取声明的样本/队列数）
#   qc             QC 收集器
# 返回：
#   临床 sample_id 向量（用于与表达矩阵交叉比对）
.check_clinical <- function(clinical_obj, cfg, qc) {

  object_name <- "clinical"

  if (!is.data.frame(clinical_obj)) {
    qc$add(
      check_name  = "clinical_type",
      object_name = object_name,
      level       = "error",
      detail      = sprintf("期望 data.frame，实际类型: %s", class(clinical_obj)[[1]])
    )
    return(character(0))
  }

  # ── 维度信息 ──
  qc$add(
    check_name  = "clinical_dims",
    object_name = object_name,
    level       = "info",
    detail      = sprintf("行数（样本）: %d，列数（字段）: %d",
                          nrow(clinical_obj), ncol(clinical_obj)),
    value       = sprintf("%d x %d", nrow(clinical_obj), ncol(clinical_obj))
  )

  # ── sample_id 列存在性 ──
  sid_col <- intersect(c("sample_id", "Sample_ID", "SampleID", "ID", "sample"),
                       names(clinical_obj))[[1]]
  if (length(sid_col) == 0) {
    qc$add(
      check_name  = "clinical_sample_id_col",
      object_name = object_name,
      level       = "error",
      detail      = "临床数据未找到 sample_id 列（候选: sample_id / Sample_ID / SampleID / ID / sample）"
    )
    return(character(0))
  }
  sample_ids <- as.character(clinical_obj[[sid_col]])

  # ── 样本 ID 唯一性 ──
  n_dup <- sum(duplicated(sample_ids))
  if (n_dup > 0) {
    qc$add(
      check_name  = "clinical_sample_uniqueness",
      object_name = object_name,
      level       = "error",
      detail      = sprintf("sample_id 存在 %d 个重复值", n_dup),
      value       = n_dup
    )
  } else {
    qc$add(
      check_name  = "clinical_sample_uniqueness",
      object_name = object_name,
      level       = "info",
      detail      = "sample_id 唯一性：通过"
    )
  }

  # ── 实际样本数 vs 声明样本数 ──
  n_actual   <- length(sample_ids)
  n_claimed  <- suppressWarnings(as.integer(cfg$n_samples_claimed))
  if (!is.na(n_claimed)) {
    deviation <- abs(n_actual - n_claimed) / n_claimed
    level_s   <- if (deviation > SAMPLE_COUNT_TOLERANCE) "warning" else "info"
    qc$add(
      check_name  = "clinical_sample_count",
      object_name = object_name,
      level       = level_s,
      detail      = sprintf(
        "实际样本数: %d，声明样本数: %d，偏差: %.1f%%",
        n_actual, n_claimed, deviation * 100
      ),
      value       = n_actual
    )
  } else {
    qc$add(
      check_name  = "clinical_sample_count",
      object_name = object_name,
      level       = "info",
      detail      = sprintf("实际样本数: %d（YAML 未声明 n_samples_claimed，跳过对比）",
                            n_actual),
      value       = n_actual
    )
  }

  # ── 队列数核查 ──
  if ("cohort" %in% names(clinical_obj)) {
    n_cohorts_actual  <- length(unique(clinical_obj$cohort))
    n_cohorts_claimed <- suppressWarnings(as.integer(cfg$n_cohorts_claimed))
    level_c <- if (!is.na(n_cohorts_claimed) &&
                   n_cohorts_actual != n_cohorts_claimed) "warning" else "info"
    detail_c <- if (!is.na(n_cohorts_claimed)) {
      sprintf("实际队列数: %d，声明队列数: %d", n_cohorts_actual, n_cohorts_claimed)
    } else {
      sprintf("实际队列数: %d（YAML 未声明 n_cohorts_claimed）", n_cohorts_actual)
    }
    qc$add(
      check_name  = "clinical_cohort_count",
      object_name = object_name,
      level       = level_c,
      detail      = detail_c,
      value       = n_cohorts_actual
    )
  } else {
    qc$add(
      check_name  = "clinical_cohort_col",
      object_name = object_name,
      level       = "warning",
      detail      = "临床数据中未找到 cohort 列，无法核查队列数"
    )
  }

  # ── 关键字段缺失率 ──
  present_key_fields <- intersect(CLINICAL_KEY_FIELDS, names(clinical_obj))
  absent_key_fields  <- setdiff(CLINICAL_KEY_FIELDS, names(clinical_obj))

  if (length(absent_key_fields) > 0) {
    qc$add(
      check_name  = "clinical_key_fields_presence",
      object_name = object_name,
      level       = "warning",
      detail      = sprintf(
        "以下关键字段在临床数据中完全缺失（列不存在）: %s",
        paste(absent_key_fields, collapse = ", ")
      )
    )
  }

  # 对已存在的关键字段逐一统计缺失率
  for (field in present_key_fields) {
    vals          <- clinical_obj[[field]]
    missing_rate  <- mean(is.na(vals) | vals %in% c("", "NA", "N/A"), na.rm = FALSE)
    level_f       <- if (missing_rate > MISSING_RATE_THRESHOLD) "warning" else "info"
    qc$add(
      check_name  = "clinical_field_completeness",
      object_name = object_name,
      level       = level_f,
      detail      = sprintf("字段 '%s' 缺失率: %.1f%%", field, missing_rate * 100),
      value       = sprintf("%.1f%%", missing_rate * 100)
    )
  }

  sample_ids
}


# ── 2c. 细胞浸润矩阵 QC ──────────────────────────────────────────────────────

# 检查细胞浸润矩阵质量（适用于 CIBERSORT / xCell / ssGSEA 输出等）
# 参数：
#   infil_obj    加载后的浸润矩阵（data.frame）
#   object_name  逻辑对象名
#   expr_sample_ids  表达矩阵的样本 ID（供交叉比对）
#   qc           QC 收集器
.check_infiltration <- function(infil_obj, object_name, expr_sample_ids, qc) {

  if (!is.data.frame(infil_obj)) {
    qc$add(
      check_name  = "infiltration_type",
      object_name = object_name,
      level       = "warning",
      detail      = sprintf("期望 data.frame，实际类型: %s（跳过浸润矩阵 QC）",
                            class(infil_obj)[[1]])
    )
    return(invisible(NULL))
  }

  # ── 维度信息 ──
  qc$add(
    check_name  = "infiltration_dims",
    object_name = object_name,
    level       = "info",
    detail      = sprintf("行数: %d，列数（细胞类型）: %d",
                          nrow(infil_obj), ncol(infil_obj)),
    value       = sprintf("%d x %d", nrow(infil_obj), ncol(infil_obj))
  )

  # ── 数值范围检查（比例值应在 0~1，分数值可能 > 1）──
  numeric_cols <- names(infil_obj)[sapply(infil_obj, is.numeric)]
  if (length(numeric_cols) > 0) {
    all_vals <- unlist(infil_obj[, numeric_cols, drop = FALSE])
    max_val  <- max(all_vals, na.rm = TRUE)
    min_val  <- min(all_vals, na.rm = TRUE)
    range_ok <- min_val >= 0 && max_val <= 1

    qc$add(
      check_name  = "infiltration_value_range",
      object_name = object_name,
      level       = if (min_val < 0) "error" else if (max_val > 1) "warning" else "info",
      detail      = sprintf(
        "数值范围: [%.4f, %.4f]。%s",
        min_val, max_val,
        if (min_val < 0) "存在负值，不符合比例值预期。"
        else if (max_val > 1) "最大值 > 1，可能为打分而非比例值，请确认。"
        else "范围正常（0~1）。"
      ),
      value = sprintf("[%.4f, %.4f]", min_val, max_val)
    )
  }

  # ── 与表达矩阵样本 ID 的交叉比对 ──
  if (length(expr_sample_ids) > 0) {
    sid_col <- intersect(c("sample_id", "Sample_ID", "SampleID", "ID"),
                         names(infil_obj))
    if (length(sid_col) > 0) {
      infil_ids <- as.character(infil_obj[[sid_col[[1]]]])
      only_in_infil <- setdiff(infil_ids,   expr_sample_ids)
      only_in_expr  <- setdiff(expr_sample_ids, infil_ids)
      level_x <- if (length(only_in_expr) > 0.1 * length(expr_sample_ids)) "warning" else "info"
      qc$add(
        check_name  = "infiltration_sample_overlap",
        object_name = object_name,
        level       = level_x,
        detail      = sprintf(
          "浸润矩阵: %d 样本；表达矩阵: %d 样本；仅在浸润矩阵: %d；仅在表达矩阵: %d",
          length(infil_ids), length(expr_sample_ids),
          length(only_in_infil), length(only_in_expr)
        )
      )
    }
  }

  invisible(NULL)
}

# ==============================================================================
# Section 3：样本 ID 交叉比对 (Cross-Object ID Check)
# ==============================================================================

# 比对表达矩阵样本 ID 与临床 sample_id 的一致性
# 参数：
#   expr_sample_ids    表达矩阵数值列名（视为样本 ID）
#   clinical_sample_ids  临床 sample_id 列
#   qc                 QC 收集器
.check_id_consistency <- function(expr_sample_ids, clinical_sample_ids, qc) {

  if (length(expr_sample_ids) == 0 || length(clinical_sample_ids) == 0) {
    qc$add(
      check_name  = "cross_id_consistency",
      object_name = "expression × clinical",
      level       = "warning",
      detail      = "样本 ID 交叉比对跳过：expression 或 clinical 样本 ID 为空"
    )
    return(invisible(NULL))
  }

  only_in_expr     <- setdiff(expr_sample_ids,     clinical_sample_ids)
  only_in_clinical <- setdiff(clinical_sample_ids, expr_sample_ids)
  n_overlap        <- length(intersect(expr_sample_ids, clinical_sample_ids))

  # 若交集为 0，极大概率是 ID 格式不一致，硬性阻断
  if (n_overlap == 0) {
    qc$add(
      check_name  = "cross_id_consistency",
      object_name = "expression × clinical",
      level       = "error",
      detail      = sprintf(
        "表达矩阵与临床数据样本 ID 无任何交集（expr: %d，clinical: %d）。\n",
        "可能原因：ID 格式不一致（如 'GSM123' vs '123'）。请手动排查。",
        length(expr_sample_ids), length(clinical_sample_ids)
      ),
      value = "overlap = 0"
    )
    return(invisible(NULL))
  }

  level_x <- if (length(only_in_expr) > 0 || length(only_in_clinical) > 0) "warning" else "info"

  qc$add(
    check_name  = "cross_id_consistency",
    object_name = "expression × clinical",
    level       = level_x,
    detail      = sprintf(
      "重叠样本: %d；仅在表达矩阵: %d；仅在临床: %d",
      n_overlap, length(only_in_expr), length(only_in_clinical)
    ),
    value       = sprintf("overlap = %d", n_overlap)
  )

  invisible(NULL)
}

# ==============================================================================
# Section 4：QC 报告写出与注册表回写 (Output)
# ==============================================================================

# 将 QC 记录写出到 CSV，并在控制台打印摘要
# 参数：
#   qc_df       QC 收集器导出的 tibble
#   dataset_id  用于文件命名
#   qc_root     输出根目录
# 返回：
#   报告文件完整路径
.write_qc_report <- function(qc_df, dataset_id, qc_root) {

  out_dir <- file.path(qc_root, QC_SUBDIR)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message(sprintf("[Info] 创建 QC 报告目录: %s", out_dir))
  }

  report_path <- file.path(out_dir, sprintf("%s_intake_qc.csv", dataset_id))
  readr::write_csv(qc_df, report_path)
  message(sprintf("[Info] QC 报告已写出: %s", report_path))

  # 控制台摘要（仅打印 warning 和 error 级别）
  notable <- dplyr::filter(qc_df, level %in% c("warning", "error"))
  if (nrow(notable) > 0) {
    message(sprintf("\n[QC 摘要] %s — 发现 %d 个需关注项:", dataset_id, nrow(notable)))
    purrr::walk2(notable$level, notable$detail, function(lvl, det) {
      prefix <- if (lvl == "error") "  ✗ [ERROR]" else "  △ [WARN]"
      message(prefix, " ", det)
    })
  }

  report_path
}

# ==============================================================================
# Section 5：核心 QC 函数 (Core QC)
# ==============================================================================

# 对单个 vendor_curated 数据集执行完整 Intake QC
#
# 参数：
#   dataset_id_or_path  dataset_id 字符串，或 YAML 配置文件路径
#   config_dir          YAML 目录，默认 CONFIG_ROOT
#   registry_path       注册表路径，默认 REGISTRY_PATH
#
# 返回：
#   list，含 status（pass/warning/fail）和 qc_report（tibble）（invisible）
#
# 副作用：
#   1. 写出 <hpc_qc_root>/intake/<dataset_id>_intake_qc.csv
#   2. 回写 intake_qc_status 和 last_update 到注册表
run_vendor_qc <- function(dataset_id_or_path,
                          config_dir    = CONFIG_ROOT,
                          registry_path = REGISTRY_PATH) {

  # ── 步骤 1：解析入参 ──
  if (file.exists(dataset_id_or_path)) {
    config_path <- normalizePath(dataset_id_or_path, mustWork = TRUE)
    dataset_id  <- tools::file_path_sans_ext(basename(config_path))
  } else {
    dataset_id  <- dataset_id_or_path
    config_path <- .resolve_config_path(dataset_id, config_dir)
  }

  message(sprintf("\n[QC] ── 开始 Intake QC: %s ──", dataset_id))

  # ── 步骤 2：读取 YAML ──
  cfg <- tryCatch(
    yaml::read_yaml(config_path),
    error = function(e) stop("[Error] YAML 解析失败: ", conditionMessage(e))
  )

  # ── 步骤 3：仅对 vendor_curated + bulk 数据集执行此脚本 ──
  if (!identical(cfg$dataset_type, "vendor_curated")) {
    stop(sprintf(
      "[Error] 此 QC 脚本仅适用于 dataset_type = vendor_curated。\n",
      "数据集 '%s' 的类型为 '%s'，请使用对应的 QC 脚本。",
      dataset_id, cfg$dataset_type
    ))
  }

  if (!identical(cfg$modality, "bulk")) {
    stop(sprintf(
      "[Error] 此 QC 脚本仅适用于 modality = bulk。\n",
      "数据集 '%s' 的 modality 为 '%s'，请使用对应的 QC 脚本。",
      dataset_id, cfg$modality
    ))
  }

  # ── 步骤 4：初始化 QC 收集器 ──
  qc <- .make_qc_collector(dataset_id)

  # ── 步骤 5：解析 raw_root 与 qc_root ──
  raw_root <- cfg$paths$raw_root
  qc_root  <- cfg$paths$qc_root %||%
              file.path(dirname(raw_root), "qc", dataset_id)

  qc$add(
    check_name  = "path_raw_root",
    object_name = "paths",
    level       = if (dir.exists(raw_root)) "info" else "error",
    detail      = sprintf("raw_root: %s（%s）",
                          raw_root,
                          if (dir.exists(raw_root)) "目录存在" else "目录不存在！"),
    value       = raw_root
  )

  if (!dir.exists(raw_root)) {
    # raw_root 不存在，无法加载任何文件，直接终止
    qc_df     <- qc$export()
    .write_qc_report(qc_df, dataset_id, qc_root)
    safe_patch(
      target_id = dataset_id,
      fields    = list(intake_qc_status = "fail",
                       last_update      = format(Sys.Date(), "%Y-%m-%d")),
      path      = registry_path
    )
    stop(sprintf("[QC FAIL] %s — raw_root 目录不存在，QC 终止。", dataset_id))
  }

  # ── 步骤 6：逐对象加载 & 执行针对性 QC ──
  files_node      <- cfg$files
  expr_sample_ids <- character(0)
  clinical_ids    <- character(0)

  for (obj_name in names(files_node)) {
    obj_def   <- files_node[[obj_name]]
    canon_fmt <- names(obj_def)[[1]]           # 第一个格式键 = canonical
    canon_file <- obj_def[[canon_fmt]]
    full_path  <- file.path(raw_root, canon_file)

    message(sprintf("  [Load] %s (%s) → %s", obj_name, canon_fmt, canon_file))

    # 加载数据对象（失败时记录 error 并跳过该对象）
    loaded_obj <- tryCatch(
      .load_data_object(full_path, canon_fmt),
      error = function(e) {
        qc$add(
          check_name  = "object_load",
          object_name = obj_name,
          level       = "error",
          detail      = sprintf("加载失败: %s", conditionMessage(e))
        )
        NULL
      }
    )

    if (is.null(loaded_obj)) next

    qc$add(
      check_name  = "object_load",
      object_name = obj_name,
      level       = "info",
      detail      = sprintf("加载成功（格式: %s）", canon_fmt)
    )

    # 根据对象名语义分流到对应 QC 检查函数
    # 命名约定（见 _template.yaml）：
    #   symbol / counts / tpm / expression → 表达矩阵
    #   clinical / metadata / sample_info  → 临床元数据
    #   cell_infiltration / infiltration   → 细胞浸润矩阵
    obj_lower <- tolower(obj_name)

    if (grepl("symbol|counts|tpm|expression|expr", obj_lower)) {
      expr_sample_ids <- .check_expression(loaded_obj, obj_name, qc)

    } else if (grepl("clinical|metadata|sample_info|phenotype", obj_lower)) {
      clinical_ids <- .check_clinical(loaded_obj, cfg, qc)

    } else if (grepl("infiltration|immune|cibersort|xcell|ssgsea", obj_lower)) {
      .check_infiltration(loaded_obj, obj_name, expr_sample_ids, qc)

    } else {
      # 未识别类型：记录 info，不执行深度检查
      qc$add(
        check_name  = "object_type_unknown",
        object_name = obj_name,
        level       = "info",
        detail      = sprintf(
          "对象 '%s' 未匹配到已知类型（expression/clinical/infiltration），跳过深度 QC。",
          obj_name
        )
      )
    }
  }

  # ── 步骤 7：样本 ID 跨对象一致性检查 ──
  .check_id_consistency(expr_sample_ids, clinical_ids, qc)

  # ── 步骤 8：汇总 QC 状态 ──
  final_status <- qc$get_status()
  qc_df        <- qc$export()

  message(sprintf(
    "\n[QC Result] %s → intake_qc_status = %s（%d checks，%d warnings，%d errors）",
    dataset_id,
    toupper(final_status),
    nrow(qc_df),
    sum(qc_df$level == "warning"),
    sum(qc_df$level == "error")
  ))

  # ── 步骤 9：写出 QC 报告 CSV ──
  .write_qc_report(qc_df, dataset_id, qc_root)

  # ── 步骤 10：回写注册表 ──
  safe_patch(
    target_id = dataset_id,
    fields    = list(
      intake_qc_status = final_status,
      hpc_qc_root      = qc_root,
      last_update      = format(Sys.Date(), "%Y-%m-%d")
    ),
    path = registry_path
  )

  # fail 时抛出可捕获的错误（批量模式下由 tryCatch 拦截）
  if (final_status == "fail") {
    stop(sprintf(
    "[QC FAIL] %s Intake QC 未通过，已写入 intake_qc_status = fail。\n详情: %s",
    dataset_id,
    detail_msg
    )) 
  }

  invisible(list(status = final_status, qc_report = qc_df))
}

# ==============================================================================
# Section 6：批量 QC 函数 (Batch QC)
# ==============================================================================

# 批量对所有 vendor_curated 数据集执行 Intake QC
#
# 参数：
#   config_dir      YAML 目录，默认 CONFIG_ROOT
#   only_pending    是否只处理 intake_qc_status = pending 的数据集（默认 TRUE）
#   stop_on_error   单个失败是否终止全局（默认 FALSE）
run_all_vendor_qc <- function(config_dir    = CONFIG_ROOT,
                              only_pending  = TRUE,
                              stop_on_error = FALSE) {

  yaml_files <- list.files(config_dir, pattern = "\\.yaml$", full.names = TRUE)
  yaml_files <- yaml_files[!grepl("/_[^/]+\\.yaml$", yaml_files)]

  # 若 only_pending = TRUE，从注册表过滤出待处理的数据集
  if (only_pending) {
    registry_df  <- read_registry(REGISTRY_PATH)
    pending_ids  <- registry_df %>%
      dplyr::filter(
        dataset_type      == "vendor_curated",
        modality          == "bulk",
        intake_qc_status  == "pending"
      ) %>%
      dplyr::pull(dataset_id)

    yaml_files <- yaml_files[
      tools::file_path_sans_ext(basename(yaml_files)) %in% pending_ids
    ]

    message(sprintf(
      "[Batch QC] 从注册表筛选出 %d 个 pending 数据集。",
      length(yaml_files)
    ))
  }

  if (length(yaml_files) == 0) {
    message("[Batch QC] 无需处理的数据集（已全部完成 QC 或目录为空）。")
    return(invisible(list(success = character(0), failed = character(0))))
  }

  success <- character(0)
  failed  <- character(0)

  for (f in yaml_files) {
    dataset_id <- tools::file_path_sans_ext(basename(f))
    result <- tryCatch(
      {
        run_vendor_qc(f)
        "ok"
      },
      error = function(e) {
        message(sprintf("[SKIP] %s QC 失败: %s", dataset_id, conditionMessage(e)))
        "error"
      }
    )

    if (result == "ok") {
      success <- c(success, dataset_id)
    } else {
      failed <- c(failed, dataset_id)
      if (stop_on_error) {
        stop(sprintf("[Fatal] 批量 QC 在 %s 处终止。", dataset_id))
      }
    }
  }

  message(sprintf(
    "\n[Batch QC Done] 成功: %d  失败: %d / 共: %d",
    length(success), length(failed), length(yaml_files)
  ))

  if (length(failed) > 0) {
    message("[Failed]:\n  ", paste(failed, collapse = "\n  "))
  }

  invisible(list(success = success, failed = failed))
}

# ==============================================================================
# Section 7：执行入口 (Execution Entry)
# ==============================================================================

if (interactive()) {
  # ── 交互式模式（RStudio Web，推荐）──────────────────────────────────────────
  # 单数据集 QC：
  # run_vendor_qc("EXAMPLE_Lung_Bulk_PublicCohort_20260101_v1")
  #
  # 批量 QC（仅处理 pending 状态）：
  # run_all_vendor_qc()
  #
  # 批量 QC（强制重跑所有，无论当前状态）：
  # run_all_vendor_qc(only_pending = FALSE)

  message(paste(
    "[Ready] 01_vendor_curated_qc.R 已加载。",
    "可用函数：",
    "  run_vendor_qc('<dataset_id 或 yaml 路径>')  ← 单数据集 QC",
    "  run_all_vendor_qc()                         ← 批量 QC（仅 pending）",
    "  run_all_vendor_qc(only_pending = FALSE)     ← 批量 QC（全部重跑）",
    sep = "\n"
  ))

} else if (sys.nframe() == 0) {
  # ── 命令行模式 ────────────────────────────────────────────────────────────
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) == 0) {
    message("[CLI] 未指定数据集，批量 QC 所有 pending 数据集...")
    run_all_vendor_qc(only_pending = TRUE, stop_on_error = FALSE)

  } else if (length(args) == 1) {
    run_vendor_qc(args[[1]])

  } else {
    stop(paste(
      "用法:",
      "  Rscript workflow/bulk/qc/01_vendor_curated_qc.R                      # 批量 QC",
      "  Rscript workflow/bulk/qc/01_vendor_curated_qc.R <dataset_id 或 path> # 单数据集",
      sep = "\n"
    ))
  }
}
