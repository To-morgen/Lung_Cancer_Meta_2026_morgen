#!/usr/bin/env Rscript
# ============================================================================
# run_normalize_pipeline.R — Phase 03 Orchestrator
# Calls: 01_merge → 02_cc_score → 03_sctransform → 04_pca → 05_cc_assess
# ============================================================================

library(here)

steps <- c(
  "01_merge_samples.R",
  "02_cell_cycle_score.R",
  "03_sctransform.R",
  "04_pca.R",
  "05_cc_assessment.R"
)

cat("\n")
cat("==============================================================\n")
cat("   Phase 03: Normalization Pipeline                            \n")
cat(sprintf("   DS_CONFIG = %s\n", Sys.getenv("DS_CONFIG", "NOT SET")))
cat(sprintf("   DS_PREFIX = %s\n", Sys.getenv("DS_PREFIX", "NOT SET")))
cat(sprintf("   Started:  %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("==============================================================\n\n")

timings <- list()
for (step in steps) {
  script <- here("workflow", "scrna", "03_normalize", step)
  cat(sprintf("\n>>> [%s] %s\n", format(Sys.time(), "%H:%M:%S"), step))
  t0 <- Sys.time()
  tryCatch({
    source(script, local = new.env())
    elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
    timings[[step]] <- list(status = "✅", time = elapsed)
    cat(sprintf(">>> ✅ %s (%.1f min)\n", step, elapsed))
  }, error = function(e) {
    elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
    timings[[step]] <<- list(status = "❌", time = elapsed)
    cat(sprintf(">>> ❌ %s FAILED (%.1f min): %s\n", step, elapsed, e$message))
  })
  gc()
}

cat("\n=== Phase 03 Summary ===\n")
for (s in names(timings)) {
  t <- timings[[s]]
  cat(sprintf("  %s %-30s %5.1f min\n", t$status, s, t$time))
}
cat(sprintf("  Finished: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

n_fail <- sum(sapply(timings, function(x) x$status == "❌"))
if (n_fail > 0) { cat(sprintf("⚠️  %d step(s) failed.\n", n_fail)); quit(status = 1) }
cat("✅ Phase 03 complete.\n")
