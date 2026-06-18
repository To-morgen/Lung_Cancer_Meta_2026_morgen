#!/bin/bash
# ============================================================================
# cellranger_utils.sh — Shared utility functions for Cell Ranger workflows
# Project: Lung_Cancer_Meta_2026
# ============================================================================

# --- Load project environment ---
load_project_env() {
    local env_file=""
    if [ -n "${PROJECT_ROOT:-}" ] && [ -f "${PROJECT_ROOT}/.env.sh" ]; then
        env_file="${PROJECT_ROOT}/.env.sh"
    else
        local script_dir
        script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
        local repo_root
        repo_root="$(cd "${script_dir}/../../.." && pwd)"
        if [ -f "${repo_root}/.env.sh" ]; then
            env_file="${repo_root}/.env.sh"
        else
            echo "[ERROR] .env.sh not found. Copy .env.sh.example to .env.sh and edit paths." >&2
            exit 1
        fi
    fi
    if [ ! -f "${env_file}" ]; then
        echo "[ERROR] .env.sh not found: ${env_file}" >&2
        exit 1
    fi
    source "${env_file}"
    echo "[env] Loaded project environment from ${env_file}"
}

# --- Validate Cell Ranger installation ---
validate_cellranger() {
    if [ ! -x "${CELLRANGER_PATH}/cellranger" ]; then
        echo "[ERROR] cellranger not found at: ${CELLRANGER_PATH}/cellranger" >&2
        exit 1
    fi
    local version
    version=$("${CELLRANGER_PATH}/cellranger" --version 2>&1 | grep -oP 'cellranger-\K[0-9.]+')
    echo "[info] Cell Ranger version: ${version}"
}

# --- Validate reference genome ---
validate_reference() {
    if [ ! -d "${CELLRANGER_REF}" ]; then
        echo "[ERROR] Reference not found: ${CELLRANGER_REF}" >&2
        exit 1
    fi
    if [ ! -f "${CELLRANGER_REF}/reference.json" ]; then
        echo "[ERROR] Invalid reference directory (missing reference.json)" >&2
        exit 1
    fi
    echo "[info] Reference: ${CELLRANGER_REF}"
}

# --- Get sample prefixes from FASTQ directory ---
# Usage: get_sample_prefixes <fastq_dir>
# Returns comma-separated list of unique sample prefixes
get_sample_prefixes() {
    local fastq_dir="$1"
    if [ ! -d "${fastq_dir}" ]; then
        echo "[ERROR] FASTQ directory not found: ${fastq_dir}" >&2
        return 1
    fi
    # Extract unique sample names from FASTQ filenames
    # Pattern: {SampleName}_S{N}_L{NNN}_R{1|2}_001.fastq.gz
    local prefixes
    prefixes=$(ls "${fastq_dir}"/*_R1_001.fastq.gz 2>/dev/null \
        | xargs -I{} basename {} \
        | sed 's/_S[0-9]*_L[0-9]*_R1_001\.fastq\.gz$//' \
        | sort -u \
        | tr '\n' ',' \
        | sed 's/,$//')
    
    if [ -z "${prefixes}" ]; then
        echo "[ERROR] No FASTQ files found in: ${fastq_dir}" >&2
        return 1
    fi
    echo "${prefixes}"
}

# --- Run cellranger count for one sample ---
# Usage: run_cellranger_count <sample_id> <fastq_dir> <output_dir> [cores] [mem_gb]
run_cellranger_count() {
    local sample_id="$1"
    local fastq_dir="$2"
    local output_dir="$3"
    local cores="${4:-32}"
    local mem_gb="${5:-128}"
    local create_bam="${6:-true}"

    # Get sample prefixes
    local prefixes
    prefixes=$(get_sample_prefixes "${fastq_dir}")
    if [ $? -ne 0 ]; then
        echo "[ERROR] Failed to get sample prefixes for ${sample_id}" >&2
        return 1
    fi

    echo "[info] Sample:    ${sample_id}"
    echo "[info] FASTQ dir: ${fastq_dir}"
    echo "[info] Prefixes:  ${prefixes}"
    echo "[info] Output:    ${output_dir}/${sample_id}/"
    echo "[info] Cores:     ${cores}, Mem: ${mem_gb}GB"

    cd "${output_dir}" || return 1

    "${CELLRANGER_PATH}/cellranger" count \
        --id="${sample_id}" \
        --fastqs="${fastq_dir}" \
        --sample="${prefixes}" \
        --transcriptome="${CELLRANGER_REF}" \
        --create-bam="${create_bam}" \
        --chemistry=auto \
        --include-introns=true \
        --localcores="${cores}" \
        --localmem="${mem_gb}"

    return $?
}

# --- Check if a sample has completed ---
# Usage: is_sample_complete <output_dir> <sample_id>
is_sample_complete() {
    local output_dir="$1"
    local sample_id="$2"
    [ -f "${output_dir}/${sample_id}/outs/metrics_summary.csv" ]
}

# --- Extract QC metrics from metrics_summary.csv ---
# Usage: extract_qc_metrics <csv_path> <sample_id>
extract_qc_metrics() {
    local csv_path="$1"
    local sample_id="$2"

    if [ ! -f "${csv_path}" ]; then
        echo "${sample_id},MISSING" >&2
        return 1
    fi

    python3 -c "
import csv, sys
with open('${csv_path}') as f:
    r = list(csv.reader(f))
    h, d = r[0], r[1]
    get = lambda name: d[h.index(name)] if name in h else 'N/A'
    keys = [
        'Estimated Number of Cells',
        'Mean Reads per Cell',
        'Median Genes per Cell',
        'Median UMI Counts per Cell',
        'Sequencing Saturation',
        'Reads Mapped Confidently to Transcriptome',
        'Fraction Reads in Cells',
        'Total Genes Detected',
        'Valid Barcodes'
    ]
    vals = [get(k) for k in keys]
    print('${sample_id},' + ','.join(vals))
"
}

echo "[info] cellranger_utils.sh loaded"
