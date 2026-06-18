#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SNAKEFILE="${ROOT}/workflow/datasets/smoke_lung_2grp/Snakefile"
SNAKEMAKE_BIN="${SNAKEMAKE_BIN:-snakemake}"

cd "${ROOT}"

Rscript examples/toy_lung_2grp/generate_toy_inputs.R
"${SNAKEMAKE_BIN}" -s "${SNAKEFILE}" --list-target-rules >/tmp/lungmeta_demo_targets.txt
"${SNAKEMAKE_BIN}" -s "${SNAKEFILE}" -n --cores 1
"${SNAKEMAKE_BIN}" -s "${SNAKEFILE}" --cores 1

echo "Demo smoke test passed."
