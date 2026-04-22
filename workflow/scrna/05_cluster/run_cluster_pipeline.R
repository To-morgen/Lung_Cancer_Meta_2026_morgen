#!/usr/bin/env Rscript
# ============================================================================
# run_cluster_pipeline.R — Phase 05 Orchestrator
# Calls: 01_cluster_umap → 02_cluster_qc
# ============================================================================

library(here)

steps <- c(
  "01_cluster_umap.R",
  "02_cluster_qc.R"
)

cat("\n")
cat("==============================================================\n")
cat("   Phase 05: Clustering Pipeline                               \n")
cat(sprintf("   DS_CONFIG = %s\n", Sys.getenv("DS_CONFIG", "NOT SET")))
cat(sprintf("   DS_PREFIX = %s\n", Sys.getenv("DS_PREFIX", "NOT SET")))
cat(sprintf("   Started:  %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("==============================================================\n\n")

timings <- list()
for (step in steps) {
  script <- here("workflow", "scrna", "05_cluster", step)
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

cat("\n=== Phase 05 Summary ===\n")
for (s in names(timings)) {
  t <- timings[[s]]
  cat(sprintf("  %s %-30s %5.1f min\n", t$status, s, t$time))
}
cat(sprintf("  Finished: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

n_fail <- sum(sapply(timings, function(x) x$status == "❌"))
if (n_fail > 0) { cat(sprintf("⚠️  %d step(s) failed.\n", n_fail)); quit(status = 1) }
cat("✅ Phase 05 complete.\n")
