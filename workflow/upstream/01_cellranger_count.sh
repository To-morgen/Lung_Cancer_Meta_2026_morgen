#!/bin/bash
# ============================================================================
# 01_cellranger_count.sh — Run Cell Ranger count for scRNA-seq samples
# Project: Lung_Cancer_Meta_2026
# Module:  upstream / alignment
#
# Usage:
#   bash workflow/upstream/01_cellranger_count.sh \
#       --dataset INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1 \
#       --samples "A1_1 A1_2 FL_1 FL_2 mc_1 mc_2" \
#       --cores 32 \
#       --mem 128 \
#       --bam true \
#       --output-subdir full
#
# Or run all samples auto-detected from fastq/ directory:
#   bash workflow/upstream/01_cellranger_count.sh \
#       --dataset INT_Novogene_LLC_scRNA_3grp6samp_20260330_v1 \
#       --auto
#
# Prerequisites:
#   - .env.sh configured with DS_RAW, CELLRANGER_PATH, CELLRANGER_REF
#   - FASTQ files organized as: {DS_RAW}/fastq/{sample_id}/
#
# Outputs:
#   - {DS_RAW}/cellranger_out/{output-subdir}/{sample_id}/outs/
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/functions/cellranger_utils.sh"

# ---- Defaults ----
CORES=32
MEM_GB=128
CREATE_BAM="true"
OUTPUT_SUBDIR="full"
AUTO_DETECT=false
SAMPLES=""
DATASET=""

# ---- Parse arguments ----
while [[ $# -gt 0 ]]; do
    case $1 in
        --dataset)    DATASET="$2";       shift 2 ;;
        --samples)    SAMPLES="$2";       shift 2 ;;
        --cores)      CORES="$2";         shift 2 ;;
        --mem)        MEM_GB="$2";        shift 2 ;;
        --bam)        CREATE_BAM="$2";    shift 2 ;;
        --output-subdir) OUTPUT_SUBDIR="$2"; shift 2 ;;
        --auto)       AUTO_DETECT=true;   shift ;;
        -h|--help)
            head -20 "$0" | grep '^#' | sed 's/^# *//'
            exit 0
            ;;
        *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# ---- Validate ----
load_project_env
validate_cellranger
validate_reference

if [ -z "${DATASET}" ]; then
    echo "[ERROR] --dataset is required" >&2
    exit 1
fi

FASTQ_BASE="${DS_RAW}/fastq"
OUTPUT_DIR="${DS_RAW}/cellranger_out/${OUTPUT_SUBDIR}"
LOG_DIR="${DS_RAW}/logs"
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

# ---- Auto-detect samples if requested ----
if [ "${AUTO_DETECT}" = true ]; then
    SAMPLES=$(ls -d "${FASTQ_BASE}"/*/ 2>/dev/null | xargs -I{} basename {} | tr '\n' ' ')
    echo "[info] Auto-detected samples: ${SAMPLES}"
fi

if [ -z "${SAMPLES}" ]; then
    echo "[ERROR] No samples specified. Use --samples or --auto" >&2
    exit 1
fi

# ---- Run ----
echo ""
echo "=========================================="
echo "Cell Ranger Count Pipeline"
echo "Dataset:  ${DATASET}"
echo "Samples:  ${SAMPLES}"
echo "Output:   ${OUTPUT_DIR}"
echo "Cores:    ${CORES}, Mem: ${MEM_GB}GB"
echo "BAM:      ${CREATE_BAM}"
echo "Started:  $(date)"
echo "=========================================="

TOTAL=0
PASSED=0
FAILED=0
SKIPPED=0

for SAMPLE in ${SAMPLES}; do
    TOTAL=$((TOTAL + 1))
    echo ""
    echo "============================================"
    echo ">>> [$(date)] Sample ${TOTAL}: ${SAMPLE}"
    echo "============================================"

    # Skip if already completed
    if is_sample_complete "${OUTPUT_DIR}" "${SAMPLE}"; then
        echo ">>> ${SAMPLE} already completed, skipping."
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # Clean up failed previous attempt
    rm -rf "${OUTPUT_DIR}/${SAMPLE}"

    # Run
    FASTQ_DIR="${FASTQ_BASE}/${SAMPLE}"
    if [ ! -d "${FASTQ_DIR}" ]; then
        echo "[ERROR] FASTQ directory not found: ${FASTQ_DIR}" >&2
        FAILED=$((FAILED + 1))
        continue
    fi

    run_cellranger_count \
        "${SAMPLE}" \
        "${FASTQ_DIR}" \
        "${OUTPUT_DIR}" \
        "${CORES}" \
        "${MEM_GB}" \
        "${CREATE_BAM}"

    if [ $? -eq 0 ]; then
        echo ">>> [$(date)] ✅ ${SAMPLE} completed successfully"
        PASSED=$((PASSED + 1))
    else
        echo ">>> [$(date)] ❌ ${SAMPLE} FAILED"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "=========================================="
echo "Pipeline finished: $(date)"
echo "Total: ${TOTAL} | Passed: ${PASSED} | Failed: ${FAILED} | Skipped: ${SKIPPED}"
echo "=========================================="

if [ ${FAILED} -gt 0 ]; then
    exit 1
fi
