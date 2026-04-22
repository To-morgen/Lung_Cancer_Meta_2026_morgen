#!/usr/bin/env Rscript
# Test if qread is actually working in the script context

suppressPackageStartupMessages({
  library(here)
  library(qs)
})

source(here("workflow", "scrna", "functions", "io_scrna.R"))

Sys.setenv(DS_PREFIX = "gse253718")

input_obj_qs  <- scrna_base("05_cluster", "objects", "seurat_clustered.qs")
input_obj_rds <- scrna_base("05_cluster", "objects", "seurat_clustered.rds")
input_obj <- if (file.exists(input_obj_qs)) input_obj_qs else input_obj_rds

cat("Selected file:", input_obj, "\n")
cat("File size:", round(file.size(input_obj) / 1e9, 1), "GB\n")
cat("File format:", tools::file_ext(input_obj), "\n\n")

cat("Starting load...\n")
t0 <- Sys.time()

if (tools::file_ext(input_obj) == "qs") {
  sobj <- qread(input_obj)
} else {
  sobj <- readRDS(input_obj)
}

t1 <- Sys.time()
cat("Load time:", round(as.numeric(difftime(t1, t0, units = "secs")), 1), "seconds\n")
cat("Object loaded:", ncol(sobj), "cells\n")
