#!/usr/bin/env Rscript
# ============================================================================
# run_annotation_pipeline.R — Phase 6: Annotation
#
# 01. FindAllMarkers (automated)
# 02. SingleR auto annotation (automated)
# ------- human reviews plots + fills CSV -------
# 03. Manual annotation from CSV (semi-automated)
# ============================================================================

library(here)

args <- commandArgs(trailingOnly = TRUE)

# Default: run steps 01 + 02 (automated part)
# To apply manual annotation: Rscript run_annotation_pipeline.R --manual
run_manual <- "--manual" %in% args

if (run_manual) {
  steps <- c("03_manual_annotate.R")
  cat("\n>>> Running manual annotation step only\n\n")
} else {
  steps <- c("01_find_markers.R", "02_singler_auto.R")
  cat("\n>>> Running automated annotation (markers + SingleR)\n")
  cat(">>> After reviewing results, run:\n")
  cat(">>>   Rscript workflow/scrna/06_annotate/run_annotation_pipeline.R --manual\n\n")
}

timings <- list()

for (step in steps) {
  script_path <- here("workflow", "scrna", "06_annotate", step)
  cat(sprintf("\n%s\n>>> [%s] %s\n%s\n\n",
              strrep("=", 60), format(Sys.time(), "%H:%M:%S"), step, strrep("=", 60)))

  t0 <- Sys.time()
  tryCatch({
    source(script_path, local = new.env())
    elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
    timings[[step]] <- list(status = "OK", time = elapsed)
    cat(sprintf("\n>>> Done: %s (%.1f min)\n", step, elapsed))
  }, error = function(e) {
    elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
    timings[[step]] <<- list(status = "FAIL", time = elapsed)
    cat(sprintf("\n>>> FAILED: %s (%.1f min): %s\n", step, elapsed, e$message))
  })
  gc()
}

cat("\n=== Phase 6 Summary ===\n")
for (s in names(timings)) {
  t <- timings[[s]]
  icon <- ifelse(t$status == "OK", "OK", "FAIL")
  cat(sprintf("  [%4s] %-30s %5.1f min\n", icon, s, t$time))
}
cat("\n")

if (!run_manual) {
  cat("NEXT STEPS:\n")
  cat("  1. Review plots:    results/scrna/06_annotate/plots/\n")
  cat("  2. Review markers:  results/scrna/06_annotate/markers/top10_markers_per_cluster.csv\n")
  cat("  3. Review SingleR:  results/scrna/06_annotate/reports/singler_cluster_annotation.csv\n")
  cat("  4. Fill mapping:    cp configs/annotation/celltype_mapping_template.csv configs/annotation/celltype_mapping.csv\n")
  cat("                      Edit celltype_mapping.csv with your annotations\n")
  cat("  5. Apply:           Rscript workflow/scrna/06_annotate/run_annotation_pipeline.R --manual\n")
}
