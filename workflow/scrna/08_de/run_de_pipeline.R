#!/usr/bin/env Rscript
# ============================================================================
# run_de_pipeline.R — Phase 08 DE Pipeline Orchestrator
# ============================================================================

library(here)

steps <- c(
  "01_pseudobulk_de.R",
  "02_enrichment.R",
  "03_de_visualization.R"
)

cat("\\n")
cat("==============================================================\\n")
cat("   Phase 08: Differential Expression Pipeline                  \\n")
cat(sprintf("   Started: %s\\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("==============================================================\\n\\n")

timings <- list()
for (step in steps) {
  script <- here("workflow", "scrna", "08_de", step)
  cat(sprintf("\\n>>> [%s] %s\\n", format(Sys.time(), "%H:%M:%S"), step))
  t0 <- Sys.time()
  tryCatch({
    source(script, local = new.env())
    elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
    timings[[step]] <- list(status = "✅", time = elapsed)
    cat(sprintf(">>> ✅ %s (%.1f min)\\n", step, elapsed))
  }, error = function(e) {
    elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
    timings[[step]] <<- list(status = "❌", time = elapsed)
    cat(sprintf(">>> ❌ %s FAILED (%.1f min): %s\\n", step, elapsed, e$message))
  })
  gc()
}

cat("\\n=== Summary ===\\n")
for (s in names(timings)) {
  t <- timings[[s]]
  cat(sprintf("  %s %-30s %5.1f min\\n", t$status, s, t$time))
}
