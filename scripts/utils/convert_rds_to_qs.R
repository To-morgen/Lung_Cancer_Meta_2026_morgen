#!/usr/bin/env Rscript
# Convert seurat_clustered.rds to qs format for faster loading

suppressPackageStartupMessages({
  library(qs)
})

input_rds <- "results/scrna/gse253718/05_cluster/objects/seurat_clustered.rds"
output_qs <- "results/scrna/gse253718/05_cluster/objects/seurat_clustered.qs"

cat("Converting RDS to QS format...\n")
cat(sprintf("  Input:  %s (%.1f GB)\n", input_rds, file.size(input_rds) / 1e9))

t0 <- Sys.time()
obj <- readRDS(input_rds)
t1 <- Sys.time()
cat(sprintf("  readRDS: %.1f seconds\n", as.numeric(difftime(t1, t0, units = "secs"))))

t2 <- Sys.time()
qsave(obj, output_qs, preset = "fast")
t3 <- Sys.time()
cat(sprintf("  qsave:   %.1f seconds\n", as.numeric(difftime(t3, t2, units = "secs"))))

cat(sprintf("  Output:  %s (%.1f GB)\n", output_qs, file.size(output_qs) / 1e9))
cat(sprintf("  Size reduction: %.1f%%\n", (1 - file.size(output_qs) / file.size(input_rds)) * 100))

# Test read speed
t4 <- Sys.time()
obj_test <- qread(output_qs)
t5 <- Sys.time()
cat(sprintf("  qread:   %.1f seconds (%.1fx faster)\n", 
            as.numeric(difftime(t5, t4, units = "secs")),
            as.numeric(difftime(t1, t0, units = "secs")) / as.numeric(difftime(t5, t4, units = "secs"))))

cat("✅ Conversion complete\n")
