#!/usr/bin/env Rscript
# ============================================================================
# run_qc_pipeline.R — Phase 2 QC Pipeline Orchestrator
#
# Usage:
#   Rscript workflow/scrna/02_qc/run_qc_pipeline.R
#
# Env vars:
#   DS_CONFIG   — path to dataset YAML
#   DS_PREFIX   — dataset output prefix
#   SKIP_SOUPX  — "true" to skip SoupX (auto-detected from config if not set)
# ============================================================================

library(here)
library(yaml)

# ── Detect whether to skip SoupX ──
skip_soupx <- FALSE

# Priority 1: explicit env var
if (tolower(Sys.getenv("SKIP_SOUPX", "")) == "true") {
  skip_soupx <- TRUE
  cat("[config] SKIP_SOUPX=true (env var)\n")
}

# Priority 2: dataset config qc_overrides
if (!skip_soupx) {
  ds_config_path <- Sys.getenv("DS_CONFIG", "")
  if (ds_config_path != "" && file.exists(ds_config_path)) {
    ds_cfg <- yaml::read_yaml(ds_config_path)
    if (isTRUE(ds_cfg$qc_overrides$skip_soupx)) {
      skip_soupx <- TRUE
      cat("[config] skip_soupx=true (from dataset YAML qc_overrides)\n")
    }
  }
}

# ── Build step list ──
all_steps <- c(
  "00_soupx_ambient.R",
  "01_create_seurat.R",
  "02_qc_filter.R",
  "03_doublet_removal.R",
  "04_qc_visualization.R"
)

if (skip_soupx) {
  steps <- all_steps[all_steps != "00_soupx_ambient.R"]
  cat("[pipeline] SoupX skipped — starting from 01_create_seurat.R\n")
} else {
  steps <- all_steps
}

cat("\n")
cat("╔══════════════════════════════════════════════════╗\n")
cat("║      Phase 2: QC Pipeline — Orchestrator        ║\n")
cat(sprintf("║  Started:  %s              ║\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("║  Steps:    %d  (SoupX: %-24s) ║\n", length(steps),
            ifelse(skip_soupx, "SKIPPED", "included")))
cat("╚══════════════════════════════════════════════════╝\n\n")

timings <- list()

for (step in steps) {
  script_path <- here("workflow", "scrna", "02_qc", step)

  cat(sprintf("\n%s\n", strrep("=", 60)))
  cat(sprintf(">>> [%s] Running: %s\n", format(Sys.time(), "%H:%M:%S"), step))
  cat(sprintf("%s\n\n", strrep("=", 60)))

  t0 <- Sys.time()

  tryCatch({
    source(script_path, local = new.env())
    elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
    timings[[step]] <- list(status = "✅", time = elapsed)
    cat(sprintf("\n>>> ✅ %s completed in %.1f min\n", step, elapsed))
  }, error = function(e) {
    elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
    timings[[step]] <<- list(status = "❌", time = elapsed)
    cat(sprintf("\n>>> ❌ %s FAILED after %.1f min: %s\n", step, elapsed, e$message))
  })

  gc()
}

# ---- Final Report ----
cat("\n\n")
cat("╔══════════════════════════════════════════════════╗\n")
cat("║            Phase 2 Pipeline Summary              ║\n")
cat("╠══════════════════════════════════════════════════╣\n")
for (s in names(timings)) {
  t <- timings[[s]]
  cat(sprintf("║  %s %-38s %5.1f min ║\n", t$status, s, t$time))
}
if (skip_soupx) {
  cat("║  ⏭️  00_soupx_ambient.R                   SKIPPED ║\n")
}
cat(sprintf("║  Finished: %s               ║\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("╚══════════════════════════════════════════════════╝\n\n")

n_fail <- sum(sapply(timings, function(x) x$status == "❌"))
if (n_fail > 0) {
  cat(sprintf("⚠️  %d step(s) failed.\n", n_fail))
  quit(status = 1)
} else {
  cat("✅ All Phase 2 steps completed successfully!\n\n")

  out_base <- scrna_base("02_qc")
  dirs <- list(
    "SoupX corrected" = file.path(out_base, "soupx"),
    "Doublet calls"   = file.path(out_base, "doublets"),
    "Filtered Seurat" = file.path(out_base, "clean"),
    "QC plots"        = file.path(out_base, "plots"),
    "Reports"         = file.path(out_base, "reports")
  )
  for (nm in names(dirs)) {
    if (dir.exists(dirs[[nm]])) {
      n <- length(list.files(dirs[[nm]], recursive = TRUE))
      cat(sprintf("  %-18s %s  (%d files)\n", nm, dirs[[nm]], n))
    }
  }
}
