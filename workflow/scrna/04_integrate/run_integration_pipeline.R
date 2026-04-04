#!/usr/bin/env Rscript
# ============================================================================
# run_integration_pipeline.R — Phase 4+5: Integration + Clustering
#
# Pipeline:
#   01. Harmony batch correction              (04_integrate/)
#   02. FindClusters + UMAP                   (05_cluster/)
#   03. Cluster QC & composition              (05_cluster/)
# ============================================================================

library(here)

steps <- list(
  list(dir = "04_integrate", script = "01_harmony.R"),
  list(dir = "05_cluster",   script = "01_cluster_umap.R"),
  list(dir = "05_cluster",   script = "02_cluster_qc.R")
)

cat("\n")
cat("╔══════════════════════════════════════════════════════╗\n")
cat("║   Phase 4+5: Integration + Clustering Pipeline      ║\n")
cat(sprintf("║   Started:  %s                ║\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("║   Steps:    %d                                        ║\n", length(steps)))
cat("╚══════════════════════════════════════════════════════╝\n\n")

timings <- list()

for (s in steps) {
  script_path <- here("workflow", "scrna", s$dir, s$script)
  label <- sprintf("%s/%s", s$dir, s$script)

  cat(sprintf("\n%s\n>>> [%s] Running: %s\n%s\n\n",
              strrep("=", 60), format(Sys.time(), "%H:%M:%S"), label, strrep("=", 60)))

  t0 <- Sys.time()
  tryCatch({
    source(script_path, local = new.env())
    elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
    timings[[label]] <- list(status = "✅", time = elapsed)
    cat(sprintf("\n>>> ✅ %s done in %.1f min\n", label, elapsed))
  }, error = function(e) {
    elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
    timings[[label]] <<- list(status = "❌", time = elapsed)
    cat(sprintf("\n>>> ❌ %s FAILED after %.1f min: %s\n", label, elapsed, e$message))
  })
  gc()
}

cat("\n\n")
cat("╔══════════════════════════════════════════════════════╗\n")
cat("║          Phase 4+5 Pipeline Summary                  ║\n")
cat("╠══════════════════════════════════════════════════════╣\n")
for (s in names(timings)) {
  t <- timings[[s]]
  cat(sprintf("║  %s %-42s %5.1f min ║\n", t$status, s, t$time))
}
cat(sprintf("║  Finished: %s                ║\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("╚══════════════════════════════════════════════════════╝\n\n")

n_fail <- sum(sapply(timings, function(x) x$status == "❌"))
if (n_fail > 0) {
  cat(sprintf("⚠️  %d step(s) failed.\n", n_fail))
  quit(status = 1)
} else {
  cat("✅ All Phase 4+5 steps completed!\n\n")
  cat("Next steps:\n")
  cat("  1. Review plots in results/scrna/05_cluster/plots/\n")
  cat("  2. Choose final clustering resolution\n")
  cat("  3. Proceed to Phase 6: Annotation\n")
}
