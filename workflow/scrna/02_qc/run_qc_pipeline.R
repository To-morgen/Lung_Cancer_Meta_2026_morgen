#!/usr/bin/env Rscript
# ============================================================================
# run_qc_pipeline.R — Phase 2 QC Pipeline Orchestrator
#
# Usage:
#   Rscript workflow/scrna/02_qc/run_qc_pipeline.R
# ============================================================================

library(here)

steps <- c(
  "01_soupx_ambient_removal.R",
  "02_scdblfinder_doublets.R",
  "03_seurat_qc_filter.R",
  "04_qc_visualization.R"
)

cat("\n")
cat("╔══════════════════════════════════════════════════╗\n")
cat("║      Phase 2: QC Pipeline — Orchestrator        ║\n")
cat(sprintf("║  Started:  %s              ║\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("║  Steps:    %d                                      ║\n", length(steps)))
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
cat(sprintf("║  Finished: %s               ║\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("╚══════════════════════════════════════════════════╝\n\n")

n_fail <- sum(sapply(timings, function(x) x$status == "❌"))
if (n_fail > 0) {
  cat(sprintf("⚠️  %d step(s) failed.\n", n_fail))
  quit(status = 1)
} else {
  cat("✅ All Phase 2 steps completed successfully!\n\n")
  
  out_base <- here("results", "scrna", "02_qc")
  dirs <- list(
    "SoupX corrected" = file.path(out_base, "soupx"),
    "Doublet calls"   = file.path(out_base, "doublets"),
    "Filtered Seurat" = file.path(out_base, "filtered"),
    "QC plots"        = file.path(out_base, "plots"),
    "Reports"         = file.path(out_base, "reports")
  )
  for (nm in names(dirs)) {
    n <- length(list.files(dirs[[nm]], recursive = TRUE))
    cat(sprintf("  %-18s %s  (%d files)\n", nm, dirs[[nm]], n))
  }
}
