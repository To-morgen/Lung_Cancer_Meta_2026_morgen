# ============================================================================
# 01_dependencies.R — Declare all project dependencies for renv detection
# 
# Purpose: renv::snapshot() scans .R files for library()/require() calls.
#          This file ensures all Phase 2+ packages are captured in renv.lock.
# ============================================================================

# --- Core analysis ---
library(Seurat)
library(SoupX)
library(scDblFinder)

# --- Bioconductor ---
library(SingleCellExperiment)
library(scran)
library(scater)
library(DropletUtils)
library(glmGamPoi)

# --- Data I/O ---
library(hdf5r)
library(rhdf5)
library(Matrix)
library(data.table)

# --- Visualization ---
library(ggplot2)
library(patchwork)

# --- Utilities ---
library(dplyr)
library(future)
library(here)
