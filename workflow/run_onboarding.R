# workflow/run_onboarding.R
#
# 职责：
#   数据集 Onboarding 全流程编排脚本（Gate 1 完整流程）。
#   按顺序调用：
#     Step 1 → 01_register_dataset.R    （注册到 dataset_registry.csv）
#     Step 2 → 02_manifest_check.R      （核查 HPC 文件完整性）
#     Step 3 → qc/01_vendor_curated_qc.R（执行 Intake QC 门控）
#
#   每步之间设置门控：前一步失败则当前数据集终止，不进入下一步。
#
# 适用类型：
#   dataset_type = vendor_curated，modality = bulk
#   其他类型数据集请使用对应的 QC 脚本（后续扩展）
#
# 三种运行模式：
#   单数据集：onboard_dataset("<dataset_id 或 yaml 路径>")
#   批量全量：onboard_all()
#   断点续跑：onboard_all(resume = TRUE)  ← 跳过已完成各步骤
#
# 版本：1.0
# 作者：Morgen

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(readr)
})

# 加载工具层（提供 read_registry / safe_patch / REGISTRY_PATH 等）
source(here("scripts",   'utils',        "utils_registry.R"))

# 加载三个步骤脚本（仅加载函数定义，不触发执行入口）
# 注：各脚本的 if(interactive()) / if(sys.nframe()==0) 块在 source 时不会执行
source(here("workflow", "intake", "01_register_dataset.R"))
source(here("workflow", "intake", "02_manifest_check.R"))
source(here("workflow", 'bulk', 'qc', '01_vendor_curated_qc.R'))

# ==============================================================================
# 常量 (Constants)
# ==============================================================================

CONFIG_ROOT <- here("configs", "datasets")

# Onboarding 报告输出路径
ONBOARDING_REPORT_PATH <- here(
  "metadata", "reports",
  sprintf("onboarding_report_%s.csv", format(Sys.Date(), "%Y%m%d"))
)

# ==============================================================================
# Section 1：内部工具 (Internal Helpers)
# ==============================================================================

# 初始化单条 onboarding 状态记录
.make_status_row <- function(dataset_id,
                             step1_register   = NA_character_,
                             step2_manifest   = NA_character_,
                             step3_intake_qc  = NA_character_,
                             overall          = NA_character_,
                             note             = NA_character_) {
  tibble::tibble(
    dataset_id       = dataset_id,
    step1_register   = step1_register,
    step2_manifest   = step2_manifest,
    step3_intake_qc  = step3_intake_qc,
    overall          = overall,
    note             = note,
    run_timestamp    = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
}

# 从注册表读取某数据集当前的步骤状态
# 用于 resume 模式判断哪些步骤已完成
.get_current_status <- function(dataset_id, registry_path = REGISTRY_PATH) {
  #' @description 查询数据集当前完成状态（用于 resume 模式判断跳过逻辑）
  #' @return list(registered, manifest, qc) 三个布尔值
  #'
  #' manifest 完成的判断逻辑（按优先级）：
  #'   1. 若注册表有 hpc_manifest_root 列且非空 → TRUE
  #'   2. 否则，检查 QC 目录下 manifest CSV 是否存在 → TRUE
  #'   3. 都不满足 → FALSE

  reg <- read_registry(registry_path)
  row <- dplyr::filter(reg, dataset_id == !!dataset_id)

  if (nrow(row) == 0) {
    return(list(registered = FALSE, manifest = FALSE, qc = FALSE))
  }

  # ── 安全读取辅助函数 ──
  .safe_get <- function(df, col) {
    if (!col %in% names(df)) return(NA_character_)
    df[[col]][[1]]
  }

  # ── manifest 完成判断 ──
  manifest_done <- FALSE

  # 方式 1：注册表中的 hpc_manifest_root 列
  manifest_val <- .safe_get(row, "hpc_manifest_root")
  if (!is.na(manifest_val) && nchar(trimws(manifest_val)) > 0) {
    manifest_done <- TRUE
  }

  # 方式 2：降级检测 — 直接检查 manifest 文件是否存在
  if (!manifest_done) {
    qc_root_val <- .safe_get(row, "hpc_qc_root")
    if (!is.na(qc_root_val)) {
      manifest_file <- file.path(
        qc_root_val, "manifest",
        paste0(dataset_id, "_manifest.csv")
      )
      manifest_done <- file.exists(manifest_file)
    }
  }

  # ── QC 完成判断 ──
  qc_val  <- .safe_get(row, "intake_qc_status")
  qc_done <- !is.na(qc_val) && qc_val %in% c("pass", "warning", "fail")

  list(
    registered = TRUE,
    manifest   = manifest_done,
    qc         = qc_done
  )
}



# 控制台分隔线（增强批量运行时的可读性）
.section_header <- function(dataset_id, step) {
  line <- strrep("─", 60)
  message(sprintf("\n%s\n[Onboarding] %s | %s\n%s", line, step, dataset_id, line))
}

# ==============================================================================
# Section 2：单数据集 Onboarding（主函数）
# ==============================================================================

# 对单个数据集执行完整的 Onboarding 流程（Step 1 → 2 → 3）
#
# 参数：
#   dataset_id_or_path  dataset_id 字符串，或 YAML 配置文件路径
#   resume              是否跳过已成功完成的步骤（默认 FALSE = 全量重跑）
#   config_dir          YAML 目录
#   registry_path       注册表路径
#
# 返回：
#   单行 tibble，记录本次 onboarding 各步骤状态（invisible）
onboard_dataset <- function(dataset_id_or_path,
                            resume        = FALSE,
                            config_dir    = CONFIG_ROOT,
                            registry_path = REGISTRY_PATH) {

  # ── 解析 dataset_id 与 config_path ──
  if (file.exists(dataset_id_or_path)) {
    config_path <- normalizePath(dataset_id_or_path, mustWork = TRUE)
    dataset_id  <- tools::file_path_sans_ext(basename(config_path))
  } else {
    dataset_id  <- dataset_id_or_path
    config_path <- file.path(config_dir, paste0(dataset_id, ".yaml"))
  }

  if (!file.exists(config_path)) {
    stop(sprintf(
      "[Error] 找不到配置文件: %s\n请先创建 YAML 配置文件。",
      config_path
    ))
  }

  # ── 初始化步骤状态 ──
  s1 <- NA_character_
  s2 <- NA_character_
  s3 <- NA_character_

  # 若 resume = TRUE，读取当前注册表状态跳过已完成步骤
  current <- if (resume) .get_current_status(dataset_id, registry_path = registry_path) else
             list(registered = FALSE, manifest = FALSE, qc = FALSE)

  # ══════════════════════════════════════════════════════════════════
  # Step 1：注册到 dataset_registry.csv
  # ════════════════════════════════════════════════
  if (resume && current$registered) {
    message(sprintf("[Skip] Step 1 Register: %s 已注册，跳过。", dataset_id))
    s1 <- "skipped(resume)"
  } else {
    .section_header(dataset_id, "Step 1 / Register")
    s1 <- tryCatch({
      register_dataset_from_yaml(config_path, registry_path = registry_path)
      "pass"
    }, error = function(e) {
      message(sprintf("[FAIL] Step 1 失败: %s", conditionMessage(e)))
      paste0("fail: ", conditionMessage(e))
    })
  }

  # Step 1 失败 → 终止
  if (startsWith(s1, "fail")) {
    return(invisible(.make_status_row(
      dataset_id,
      step1_register  = s1,
      overall         = "fail",
      note            = "Step 1 Register 失败，流程终止"
    )))
  }

  # ══════════════════════════════════════════════════════════════════
  # Step 2：Manifest 核查
  # ══════════════════════════════════════════════════════════════════
  if (resume && current$manifest) {
    message(sprintf("[Skip] Step 2 Manifest: %s 已完成，跳过。", dataset_id))
    s2 <- "skipped(resume)"
  } else {
    .section_header(dataset_id, "Step 2 / Manifest Check")
    s2 <- tryCatch({
      check_dataset_manifest(config_path, registry_path = registry_path)
      "pass"
    }, error = function(e) {
      message(sprintf("[FAIL] Step 2 失败: %s", conditionMessage(e)))
      paste0("fail: ", conditionMessage(e))
    })
  }

  # Step 2 失败 → 终止（canonical 文件缺失）
  if (startsWith(s2, "fail")) {
    return(invisible(.make_status_row(
      dataset_id,
      step1_register  = s1,
      step2_manifest  = s2,
      overall         = "fail",
      note            = "Step 2 Manifest 失败，请确认文件已完整下载"
    )))
  }

  # ══════════════════════════════════════════════════════════════════
  # Step 3：Intake QC（Gate 1）
  # ══════════════════════════════════════════════════════════════════
  if (resume && current$qc) {
    # resume 模式下，QC 已完成（pass/warning/fail 均视为"已跑"，不重跑）
    reg_row <- dplyr::filter(read_registry(registry_path), dataset_id == !!dataset_id)
    s3 <- paste0("skipped(resume):", reg_row$intake_qc_status[[1]])
    message(sprintf("[Skip] Step 3 QC: %s 已完成（%s），跳过。",
                    dataset_id, reg_row$intake_qc_status[[1]]))
  } else {
    .section_header(dataset_id, "Step 3 / Intake QC")
    s3 <- tryCatch({
      result <- run_vendor_qc(config_path, registry_path = registry_path)
      result$status   # "pass" 或 "warning"
    }, error = function(e) {
      # run_vendor_qc 在 fail 时会 stop()，此处捕获
      message(sprintf("[FAIL] Step 3 QC: %s", conditionMessage(e)))
      "fail"
    })
  }

  # ── 汇总 overall 状态 ──
  overall <- dplyr::case_when(
    startsWith(s1, "fail") | startsWith(s2, "fail") | s3 == "fail" ~ "fail",
    s3 == "warning" | grepl("warning", s3) ~ "warning",
    TRUE ~ "pass"
  )

  line <- strrep("═", 60)
  message(sprintf(
    "\n%s\n[Onboarding Done] %s\n  Step 1 Register : %s\n  Step 2 Manifest : %s\n  Step 3 QC       : %s\n  Overall         : %s\n%s",
    line, dataset_id, s1, s2, s3, toupper(overall), line
  ))

  invisible(.make_status_row(
    dataset_id,
    step1_register  = s1,
    step2_manifest  = s2,
    step3_intake_qc = s3,
    overall         = overall
  ))
}

# ==============================================================================
# Section 3：批量 Onboarding
# ==============================================================================

# 批量对 configs/datasets/ 下所有数据集执行完整 Onboarding 流程
#
# 参数：
#   config_dir     YAML 目录
#   resume         TRUE = 跳过已完成步骤（断点续跑）；FALSE = 全量重跑
#   stop_on_error  TRUE = 遇到失败立即终止全局；FALSE = 跳过并记录（推荐）
#   dataset_types  限定处理的 dataset_type（默认仅处理 vendor_curated）
#
# 返回：
#   汇总报告 tibble（同时写出到 ONBOARDING_REPORT_PATH）
onboard_all <- function(config_dir    = CONFIG_ROOT,
                        resume        = FALSE,
                        stop_on_error = FALSE,
                        dataset_types = "vendor_curated") {

  yaml_files <- list.files(config_dir, pattern = "\\.yaml$", full.names = TRUE)
  # 排除模板文件
  yaml_files <- yaml_files[!grepl("/_[^/]+\\.yaml$", yaml_files)]

  if (length(yaml_files) == 0) {
    message("[Warning] 未找到可处理的 YAML 配置文件。")
    return(invisible(tibble::tibble()))
  }

  # 若指定了 dataset_types，仅保留匹配类型的 YAML
  # （通过快速读取 YAML 顶层字段实现，不依赖注册表）
  if (!is.null(dataset_types)) {
    yaml_files <- purrr::keep(yaml_files, function(f) {
      cfg <- tryCatch(yaml::read_yaml(f), error = function(e) NULL)
      if (is.null(cfg)) return(FALSE)
      cfg$dataset_type %in% dataset_types
    })
    message(sprintf(
      "[Batch] 过滤后剩余 %d 个 %s 数据集。",
      length(yaml_files), paste(dataset_types, collapse = "/")
    ))
  }

  if (length(yaml_files) == 0) {
    message("[Warning] 过滤后无匹配数据集。")
    return(invisible(tibble::tibble()))
  }

  message(sprintf(
    "\n[Batch Onboarding] 共 %d 个数据集，resume = %s\n",
    length(yaml_files), resume
  ))

  report_rows <- purrr::map(yaml_files, function(f) {
    dataset_id <- tools::file_path_sans_ext(basename(f))
    row <- tryCatch(
      onboard_dataset(f, resume = resume),
      error = function(e) {
        message(sprintf("[SKIP] %s 意外错误: %s", dataset_id, conditionMessage(e)))
        .make_status_row(
          dataset_id,
          overall = "error",
          note    = conditionMessage(e)
        )
      }
    )

    if (!is.null(row) && row$overall == "fail" && stop_on_error) {
      stop(sprintf("[Fatal] Batch onboarding 在 %s 处终止（stop_on_error = TRUE）。",
                   dataset_id))
    }

    row
  })

  report_df <- dplyr::bind_rows(report_rows)

  # ── 写出汇总报告 ──
  report_dir <- dirname(ONBOARDING_REPORT_PATH)
  if (!dir.exists(report_dir)) {
    dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  }
  readr::write_csv(report_df, ONBOARDING_REPORT_PATH)

  # ── 控制台汇总 ──
  n_pass    <- sum(report_df$overall == "pass",    na.rm = TRUE)
  n_warning <- sum(report_df$overall == "warning", na.rm = TRUE)
  n_fail    <- sum(report_df$overall %in% c("fail", "error"), na.rm = TRUE)

  message(sprintf(
    "\n%s\n[Batch Onboarding Done]\n  ✓ pass    : %d\n  △ warning : %d\n  ✗ fail    : %d\n  总计      : %d\n  报告路径  : %s\n%s",
    strrep("═", 60),
    n_pass, n_warning, n_fail, nrow(report_df),
    ONBOARDING_REPORT_PATH,
    strrep("═", 60)
  ))

  invisible(report_df)
}

# ==============================================================================
# Section 4：快捷查询函数 (Quick Inspection)
# ==============================================================================

# 打印当前所有数据集的 Onboarding 状态概览（从注册表读取）
show_onboarding_status <- function(registry_path = REGISTRY_PATH) {
  reg <- read_registry(registry_path)

  if (nrow(reg) == 0) {
    message("[Info] 注册表为空。")
    return(invisible(NULL))
  }

  # 选择关键状态列展示（若列不存在则跳过）
  key_cols <- intersect(
    c("dataset_id", "dataset_type", "modality", "status",
      "intake_qc_status", "harmonization_status", "last_update"),
    names(reg)
  )

  summary_df <- dplyr::select(reg, dplyr::all_of(key_cols)) %>%
    dplyr::arrange(intake_qc_status, dataset_id)

  print(summary_df, n = Inf)
  invisible(summary_df)
}

# ==============================================================================
# Section 5：执行入口 (Execution Entry)
# ==============================================================================

if (interactive()) {
  # ── 交互式模式（RStudio Web，推荐）──────────────────────────────────────────
  #
  # 单数据集全量 Onboarding：
  # onboard_dataset("VEND_AIshixin_LUAD_Bulk_22Cohorts_3019_20260307_v1")
  #
  # 单数据集断点续跑（跳过已完成步骤）：
  # onboard_dataset("VEND_AIshixin_LUAD_Bulk_22Cohorts_3019_20260307_v1",
  #                 resume = TRUE)
  #
  # 批量全量 Onboarding：
  # onboard_all()
  #
  # 批量断点续跑：
  # onboard_all(resume = TRUE)
  #
  # 查看所有数据集当前状态：
  # show_onboarding_status()

  message(paste(
    "[Ready] run_onboarding.R 已加载。",
    "可用函数：",
    "  onboard_dataset('<dataset_id 或 yaml>')       ← 单数据集全量",
    "  onboard_dataset('...', resume = TRUE)         ← 单数据集断点续跑",
    "  onboard_all()                                 ← 批量全量",
    "  onboard_all(resume = TRUE)                    ← 批量断点续跑",
    "  show_onboarding_status()                      ← 查看当前状态概览",
    sep = "\n"
  ))

} else if (sys.nframe() == 0) {
  # ── 命令行模式 ────────────────────────────────────────────────────────────
  args <- commandArgs(trailingOnly = TRUE)

  # 支持可选标志：--resume / --stop-on-error
  resume_flag        <- "--resume"        %in% args
  stop_on_error_flag <- "--stop-on-error" %in% args
  pos_args           <- args[!args %in% c("--resume", "--stop-on-error")]

  if (length(pos_args) == 0) {
    # 无位置参数：批量模式
    message(sprintf(
      "[CLI] 批量 Onboarding（resume=%s，stop_on_error=%s）...",
      resume_flag, stop_on_error_flag
    ))
    onboard_all(resume = resume_flag, stop_on_error = stop_on_error_flag)

  } else if (length(pos_args) == 1) {
    # 一个位置参数：单数据集
    onboard_dataset(pos_args[[1]], resume = resume_flag)

  } else {
    stop(paste(
      "用法:",
      "  Rscript workflow/run_onboarding.R [--resume] [--stop-on-error]",
      "  Rscript workflow/run_onboarding.R <dataset_id 或 yaml> [--resume]",
      sep = "\n"
    ))
  }
}
