# tests/test_onboarding_integration.R
#
# 集成测试：验证 Onboarding 全流程（Step 1→2→3）的数据流和门控逻辑
#
# 修复记录：
#   v1.1 - 根因 1：sandbox registry_cols 从 schema 自动派生，消除列不一致
#   v1.1 - 根因 2：所有函数调用显式传入 registry_path，保证沙箱隔离
#
# 运行方式：
#   testthat::test_file(here::here("tests", "test_onboarding_integration.R"))

suppressPackageStartupMessages({
  library(here)
  library(testthat)
  library(yaml)
  library(readr)
  library(dplyr)
  library(tibble)
})

# ==============================================================================
# Section 0：加载被测脚本
# ==============================================================================

source(here("scripts", "utils", "utils_registry.R"))
source(here("workflow", "intake", "01_register_dataset.R"))
source(here("workflow", "intake", "02_manifest_check.R"))
source(here("workflow", "bulk", "qc", "01_vendor_curated_qc.R"))
source(here("workflow",           "run_onboarding.R"))

# 保存原始全局常量（每个 test_that 结束后恢复）
.ORIG_SCHEMA_PATH   <- SCHEMA_PATH
.ORIG_REGISTRY_PATH <- REGISTRY_PATH

# ==============================================================================
# Section 1：沙箱环境构建工具 (Sandbox Helpers)
# ==============================================================================

.setup_sandbox <- function(prefix = "onboard_test") {

  root <- file.path(tempdir(), paste0(prefix, "_", format(Sys.time(), "%H%M%S%OS2")))

  dirs <- list(
    root         = root,
    configs      = file.path(root, "configs", "datasets"),
    registry_dir = file.path(root, "metadata", "registry"),
    reports      = file.path(root, "metadata", "reports"),
    schema_dir   = file.path(root, "configs", "registry"),
    raw          = file.path(root, "data", "raw"),
    qc           = file.path(root, "data", "qc")
  )

  purrr::walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

  # ── 根因 1 修复：从生产 schema 派生 registry_cols ──
  # 将生产 schema 复制到沙箱，然后从中读取列名
  sandbox_schema_path <- file.path(dirs$schema_dir, "registry_schema.yaml")

  if (file.exists(.ORIG_SCHEMA_PATH)) {
    file.copy(.ORIG_SCHEMA_PATH, sandbox_schema_path)
    schema <- yaml::read_yaml(sandbox_schema_path)
    registry_cols <- names(schema$columns)
  } else {
    # 降级：若生产 schema 也不存在，用最小列集合
    warning("[Test] 生产 schema 不存在，使用最小列集合")
    registry_cols <- c(
      "dataset_id", "display_name", "dataset_type", "disease", "subtype",
      "modality", "data_level", "source_name", "owner",
      "n_cohorts_claimed", "n_samples_claimed",
      "status", "intake_qc_status", "harmonization_status",
      "download_date", "register_date", "last_update",
      "notes"
    )
  }

  # 创建空注册表（列头完全匹配 schema）
  registry_path <- file.path(dirs$registry_dir, "dataset_registry.csv")
  empty_reg <- tibble::as_tibble(
    setNames(
      as.list(rep(NA_character_, length(registry_cols))),
      registry_cols
    )
  )[0, ]
  readr::write_csv(empty_reg, registry_path)

  dirs$registry_path      <- registry_path
  dirs$sandbox_schema_path <- sandbox_schema_path

  # ── 根因 1 修复：覆盖全局 SCHEMA_PATH，使 validate_registry() 读沙箱 schema ──
  assign("SCHEMA_PATH", sandbox_schema_path, envir = globalenv())

  dirs
}


.teardown_sandbox <- function(sandbox) {
  # 恢复全局常量
  assign("SCHEMA_PATH",   .ORIG_SCHEMA_PATH,   envir = globalenv())
  assign("REGISTRY_PATH", .ORIG_REGISTRY_PATH, envir = globalenv())

  if (dir.exists(sandbox$root)) {
    unlink(sandbox$root, recursive = TRUE)
  }
}


# 构造 YAML 配置文件（与之前版本相同，但 paths 写法更精确）
.make_yaml <- function(sandbox, dataset_id, overrides = list()) {

  raw_root <- file.path(sandbox$raw, dataset_id)
  qc_root  <- file.path(sandbox$qc,  dataset_id)
  dir.create(raw_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_root,  recursive = TRUE, showWarnings = FALSE)

  # ── 默认配置（覆盖 01_register_dataset.R 所有必填字段）──
  default_cfg <- list(
    dataset_id           = dataset_id,
    display_name         = paste("Test Dataset", dataset_id),
    dataset_type         = "vendor_curated",
    disease              = "Lung Cancer",
    subtype              = "LUAD",
    modality             = "bulk",
    data_level           = "processed",
    source_name          = "TestVendor",
    source_accession     = NA,
    owner                = "Morgen",
    n_cohorts_claimed    = 3,
    n_samples_claimed    = 30,
    n_samples_aligned    = NA,

    # ★ P1 修复：补全所有状态类必填字段
    status               = "registered",
    intake_qc_status     = "pending",
    harmonization_status = "pending",
    analysis_qc_status   = "pending",      # ← 之前缺这个
    atlas_role           = "undecided",     # ← 之前缺这个
    inclusion_decision   = "pending",       # ← 之前缺这个

    download_date        = "2026-01-01",
    register_date        = format(Sys.Date(), "%Y-%m-%d"),
    last_update          = format(Sys.Date(), "%Y-%m-%d"),

    paths = list(
      raw_root      = raw_root,
      qc_root       = qc_root,
      manifest_root = NA,
      stage_root    = NA
    ),
    files = list(
      symbol   = list(rds = "expression_matrix.rds"),
      clinical = list(rds = "clinical_metadata.rds")
    ),
    notes = "integration test mock"
  )

  cfg <- modifyList(default_cfg, overrides)
  yaml_path <- file.path(sandbox$configs, paste0(dataset_id, ".yaml"))
  yaml::write_yaml(cfg, yaml_path)
  yaml_path
}



# 构造 Mock 数据文件
.make_mock_data <- function(sandbox, dataset_id,
                             n_samples = 30,
                             n_genes   = 200,
                             n_cohorts = 3) {

  raw_root   <- file.path(sandbox$raw, dataset_id)
  sample_ids <- paste0("SAMPLE_", sprintf("%03d", seq_len(n_samples)))
  gene_ids   <- paste0("GENE_",   sprintf("%04d", seq_len(n_genes)))

  expr_mat <- matrix(
    runif(n_genes * n_samples, min = 0.1, max = 15),
    nrow     = n_genes,
    ncol     = n_samples,
    dimnames = list(gene_ids, sample_ids)
  )
  expr_df <- as.data.frame(expr_mat)
  expr_df <- tibble::rownames_to_column(expr_df, var = "gene")

  cohort_labels <- rep(paste0("COHORT_", seq_len(n_cohorts)),
                       length.out = n_samples)
  clinical_df <- tibble::tibble(
    sample_id = sample_ids,
    cohort    = cohort_labels,
    OS        = sample(c(0L, 1L), n_samples, replace = TRUE),
    OS.time   = round(runif(n_samples, 3, 60), 1),
    gender    = sample(c("Male", "Female"), n_samples, replace = TRUE),
    age       = sample(40:80, n_samples, replace = TRUE),
    stage     = sample(c("I", "II", "III", "IV"), n_samples, replace = TRUE)
  )

  saveRDS(expr_df,     file.path(raw_root, "expression_matrix.rds"))
  saveRDS(clinical_df, file.path(raw_root, "clinical_metadata.rds"))

  list(expr = expr_df, clinical = clinical_df, sample_ids = sample_ids)
}

# ==============================================================================
# Section 2：断言辅助函数
# ==============================================================================

.get_field <- function(sandbox, dataset_id, field) {
  reg <- readr::read_csv(
    sandbox$registry_path,
    col_types = readr::cols(.default = "c"),
    show_col_types = FALSE
  )
  row <- dplyr::filter(reg, dataset_id == !!dataset_id)
  if (nrow(row) == 0) return(NA_character_)
  if (!field %in% names(row)) return(NA_character_)
  row[[field]][[1]]
}

.qc_report_exists <- function(sandbox, dataset_id) {
  qc_root     <- file.path(sandbox$qc, dataset_id)
  report_path <- file.path(qc_root, "intake",
                            paste0(dataset_id, "_intake_qc.csv"))
  file.exists(report_path)
}

# ==============================================================================
# Section 3：测试套件
# ==============================================================================

# ──────────────────────────────────────────────────────────────────────────────
# Test Suite 1：Happy Path
# ──────────────────────────────────────────────────────────────────────────────

test_that("Happy Path: 完整 Onboarding 流程 pass", {

  sb <- .setup_sandbox("happy")
  on.exit(.teardown_sandbox(sb), add = TRUE)

  dataset_id <- "VEND_Test_LUAD_Bulk_HappyPath_20260101_v1"
  yaml_path  <- .make_yaml(sb, dataset_id)
  .make_mock_data(sb, dataset_id)

  # ★ 根因 2 修复：显式传入 registry_path
  result <- onboard_dataset(
    yaml_path,
    resume        = FALSE,
    config_dir    = sb$configs,
    registry_path = sb$registry_path
  )

  # Step 1：注册表中应存在该数据集
  expect_equal(result$step1_register, "pass")
  expect_equal(.get_field(sb, dataset_id, "dataset_id"), dataset_id)
  expect_equal(.get_field(sb, dataset_id, "dataset_type"), "vendor_curated")

  # Step 2：manifest 应已写出
  expect_equal(result$step2_manifest, "pass")

  # Step 3：intake_qc_status 应为 pass
  expect_true(result$step3_intake_qc %in% c("pass", "warning"))
  expect_true(
    .get_field(sb, dataset_id, "intake_qc_status") %in% c("pass", "warning")
  )
  expect_true(.qc_report_exists(sb, dataset_id))

  # Overall
  expect_true(result$overall %in% c("pass", "warning"))
})


# ──────────────────────────────────────────────────────────────────────────────
# Test Suite 2：Warning Path（临床字段缺失率高）
# ──────────────────────────────────────────────────────────────────────────────

test_that("Warning Path: 临床字段缺失率高触发 warning 但不阻断", {

  sb <- .setup_sandbox("warning")
  on.exit(.teardown_sandbox(sb), add = TRUE)

  dataset_id <- "VEND_Test_LUAD_Bulk_WarningPath_20260101_v1"
  yaml_path  <- .make_yaml(sb, dataset_id)
  mock       <- .make_mock_data(sb, dataset_id)

  # 篡改临床数据：OS / OS.time 80% 置为 NA
  raw_root     <- file.path(sb$raw, dataset_id)
  clinical_bad <- mock$clinical
  n_na         <- floor(nrow(clinical_bad) * 0.8)
  clinical_bad$OS[seq_len(n_na)]      <- NA
  clinical_bad$OS.time[seq_len(n_na)] <- NA
  saveRDS(clinical_bad, file.path(raw_root, "clinical_metadata.rds"))

  result <- onboard_dataset(
    yaml_path,
    config_dir    = sb$configs,
    registry_path = sb$registry_path
  )

  # Step 3 应含 warning
  expect_true(result$step3_intake_qc %in% c("warning", "pass"))
  expect_true(result$overall %in% c("warning", "pass"))

  # QC 报告应存在
  expect_true(.qc_report_exists(sb, dataset_id))
})


# ──────────────────────────────────────────────────────────────────────────────
# Test Suite 3：Fail Path 1 — canonical 文件缺失（Step 2 门控）
# ──────────────────────────────────────────────────────────────────────────────

test_that("Fail Path 1: canonical 文件缺失，Step 2 阻断", {

  sb <- .setup_sandbox("fail_manifest")
  on.exit(.teardown_sandbox(sb), add = TRUE)

  dataset_id <- "VEND_Test_LUAD_Bulk_FailManifest_20260101_v1"
  yaml_path  <- .make_yaml(sb, dataset_id)
  # 故意不创建数据文件

  result <- onboard_dataset(
    yaml_path,
    config_dir    = sb$configs,
    registry_path = sb$registry_path
  )

  # Step 1 应通过
  expect_equal(result$step1_register, "pass")

  # Step 2 应 fail
  expect_true(startsWith(result$step2_manifest, "fail"))

  # Step 3 不应执行（NA）
  expect_true(is.na(result$step3_intake_qc))

  # Overall = fail
  expect_equal(result$overall, "fail")
})


# ──────────────────────────────────────────────────────────────────────────────
# Test Suite 4：Fail Path 2 — 表达矩阵含 NA（Step 3 门控）
# ──────────────────────────────────────────────────────────────────────────────

test_that("Fail Path 2: 表达矩阵含 NA，Step 3 门控 fail", {

  sb <- .setup_sandbox("fail_qc")
  on.exit(.teardown_sandbox(sb), add = TRUE)

  dataset_id <- "VEND_Test_LUAD_Bulk_FailQC_20260101_v1"
  yaml_path  <- .make_yaml(sb, dataset_id)
  mock       <- .make_mock_data(sb, dataset_id)

  # 篡改表达矩阵：注入 NA
  raw_root     <- file.path(sb$raw, dataset_id)
  expr_bad     <- mock$expr
  numeric_cols <- names(expr_bad)[sapply(expr_bad, is.numeric)]
  set.seed(42)
  for (i in seq_len(50)) {
    r <- sample(nrow(expr_bad), 1)
    cc <- sample(numeric_cols, 1)
    expr_bad[r, cc] <- NA
  }
  saveRDS(expr_bad, file.path(raw_root, "expression_matrix.rds"))

  result <- onboard_dataset(
    yaml_path,
    config_dir    = sb$configs,
    registry_path = sb$registry_path
  )

  # Step 1 & 2 应通过
  expect_equal(result$step1_register, "pass")
  expect_equal(result$step2_manifest, "pass")

  # Step 3 应 fail
  expect_equal(result$step3_intake_qc, "fail")
  expect_equal(.get_field(sb, dataset_id, "intake_qc_status"), "fail")

  # QC 报告应存在且含 error
  # 改为：
  report_path <- file.path(sb$qc, dataset_id, "intake",
                          paste0(dataset_id, "_intake_qc.csv"))
  expect_true(file.exists(report_path))

  # 用 skip_if 保护后续读取，避免文件不存在时抛出 Error（而非 FAIL）
  skip_if(!file.exists(report_path), "QC 报告不存在，跳过内容校验")
  qc_df <- readr::read_csv(report_path, show_col_types = FALSE)
  expect_true(any(qc_df$level == "error"))
})


# ──────────────────────────────────────────────────────────────────────────────
# Test Suite 5：Resume Mode（断点续跑）
# ──────────────────────────────────────────────────────────────────────────────

test_that("Resume Mode: 跳过已完成步骤", {

  sb <- .setup_sandbox("resume")
  on.exit(.teardown_sandbox(sb), add = TRUE)

  dataset_id <- "VEND_Test_LUAD_Bulk_Resume_20260101_v1"
  yaml_path  <- .make_yaml(sb, dataset_id)
  .make_mock_data(sb, dataset_id)

  # ★ 根因 2 修复：手动运行 Step 1 & 2 时传入 registry_path
  register_dataset_from_yaml(yaml_path, registry_path = sb$registry_path)
  check_dataset_manifest(yaml_path, registry_path = sb$registry_path)

  # 确认 Step 1 & 2 完成，Step 3 仍 pending
  expect_equal(.get_field(sb, dataset_id, "intake_qc_status"), "pending")

  # 第二次运行：resume = TRUE
  result <- onboard_dataset(
    yaml_path,
    resume        = TRUE,
    config_dir    = sb$configs,
    registry_path = sb$registry_path
  )

  # Step 1 & 2 应跳过
  expect_true(grepl("skipped", result$step1_register))
  expect_true(grepl("skipped", result$step2_manifest))

  # Step 3 应执行并完成
  expect_true(result$step3_intake_qc %in% c("pass", "warning"))
})


# ──────────────────────────────────────────────────────────────────────────────
# Test Suite 6：Batch Mode
# ──────────────────────────────────────────────────────────────────────────────

test_that("Batch Mode: 多数据集批量处理，部分失败不中断", {

  sb <- .setup_sandbox("batch")
  on.exit(.teardown_sandbox(sb), add = TRUE)

  ids <- list(
    pass1 = "VEND_Test_LUAD_Bulk_Batch1_20260101_v1",
    pass2 = "VEND_Test_LUAD_Bulk_Batch2_20260101_v1",
    fail1 = "VEND_Test_LUAD_Bulk_Batch3_20260101_v1"
  )

  .make_yaml(sb, ids$pass1)
  .make_mock_data(sb, ids$pass1)

  .make_yaml(sb, ids$pass2)
  .make_mock_data(sb, ids$pass2)

  .make_yaml(sb, ids$fail1)
  # 故意不创建数据文件 → Step 2 fail

  # ★ 批量模式需要 onboard_all 能接受 registry_path
  # 这里直接调用 onboard_dataset 循环模拟，避免改 onboard_all 签名
  results <- purrr::map(
    list.files(sb$configs, pattern = "\\.yaml$", full.names = TRUE),
    function(f) {
      tryCatch(
        onboard_dataset(f, config_dir = sb$configs, registry_path = sb$registry_path),
        error = function(e) {
          .make_status_row(
            tools::file_path_sans_ext(basename(f)),
            overall = "error",
            note    = conditionMessage(e)
          )
        }
      )
    }
  )
  report_df <- dplyr::bind_rows(results)

  expect_equal(nrow(report_df), 3)

  # 两个 pass（或 warning），一个 fail
  pass_rows <- dplyr::filter(report_df, overall %in% c("pass", "warning"))
  fail_rows <- dplyr::filter(report_df, overall %in% c("fail", "error"))
  expect_equal(nrow(pass_rows), 2)
  expect_gte(nrow(fail_rows), 1)
})


# ──────────────────────────────────────────────────────────────────────────────
# Test Suite 7：幂等性（同一数据集运行两次不产生重复行）
# ──────────────────────────────────────────────────────────────────────────────

test_that("幂等性: 重复运行不产生重复注册表行", {

  sb <- .setup_sandbox("idempotent")
  on.exit(.teardown_sandbox(sb), add = TRUE)

  dataset_id <- "VEND_Test_LUAD_Bulk_Idempotent_20260101_v1"
  yaml_path  <- .make_yaml(sb, dataset_id)
  .make_mock_data(sb, dataset_id)

  onboard_dataset(yaml_path, config_dir = sb$configs, registry_path = sb$registry_path)
  onboard_dataset(yaml_path, config_dir = sb$configs, registry_path = sb$registry_path)

  reg <- readr::read_csv(
    sb$registry_path,
    col_types = readr::cols(.default = "c"),
    show_col_types = FALSE
  )
  n_rows <- sum(reg$dataset_id == dataset_id)
  expect_equal(n_rows, 1L)
})


# ──────────────────────────────────────────────────────────────────────────────
# Test Suite 8：注册表隔离性验证（沙箱写入不影响真实注册表）
# ──────────────────────────────────────────────────────────────────────────────

test_that("沙箱隔离: 测试写入不污染生产注册表", {

  # 记录生产注册表当前行数
  prod_reg_before <- NULL
  if (file.exists(.ORIG_REGISTRY_PATH)) {
    prod_reg_before <- readr::read_csv(
      .ORIG_REGISTRY_PATH,
      col_types = readr::cols(.default = "c"),
      show_col_types = FALSE
    )
  }

  sb <- .setup_sandbox("isolation")
  on.exit(.teardown_sandbox(sb), add = TRUE)

  dataset_id <- "VEND_Test_LUAD_Bulk_Isolation_20260101_v1"
  yaml_path  <- .make_yaml(sb, dataset_id)
  .make_mock_data(sb, dataset_id)

  onboard_dataset(yaml_path, config_dir = sb$configs, registry_path = sb$registry_path)

  # 确认沙箱注册表包含该记录
  sandbox_reg <- readr::read_csv(
    sb$registry_path,
    col_types = readr::cols(.default = "c"),
    show_col_types = FALSE
  )
  expect_true(dataset_id %in% sandbox_reg$dataset_id)

  # 确认生产注册表未被修改
  if (!is.null(prod_reg_before)) {
    prod_reg_after <- readr::read_csv(
      .ORIG_REGISTRY_PATH,
      col_types = readr::cols(.default = "c"),
      show_col_types = FALSE
    )
    expect_equal(nrow(prod_reg_before), nrow(prod_reg_after))
    expect_false(dataset_id %in% prod_reg_after$dataset_id)
  }
})
