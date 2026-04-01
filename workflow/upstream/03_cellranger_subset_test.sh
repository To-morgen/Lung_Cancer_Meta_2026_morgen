#!/bin/bash
# ============================================================================
# 03_cellranger_subset_test.sh — Create FASTQ subset and run Cell Ranger test
#
# Usage:
#   bash workflow/upstream/03_cellranger_subset_test.sh \
#       --sample A1_1 \
#       --n-reads 1000000 \
#       --cores 8 \
#       --mem 32
#
# Purpose:
#   Validate Cell Ranger setup (reference, chemistry, FASTQ naming) before
#   committing to full run. Uses seqtk to subsample reads.
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/functions/cellranger_utils.sh"

SAMPLE=""
N_READS=1000000
CORES=8
MEM_GB=32
SEED=42

while [[ $# -gt 0 ]]; do
    case $1 in
        --sample)   SAMPLE="$2";  shift 2 ;;
        --n-reads)  N_READS="$2"; shift 2 ;;
        --cores)    CORES="$2";   shift 2 ;;
        --mem)      MEM_GB="$2";  shift 2 ;;
        --seed)     SEED="$2";    shift 2 ;;
        -h|--help)  head -15 "$0" | grep '^#' | sed 's/^# *//'; exit 0 ;;
        *)          echo "[ERROR] Unknown: $1" >&2; exit 1 ;;
    esac
done

load_project_env
validate_cellranger
validate_reference

if [ -z "${SAMPLE}" ]; then
    echo "[ERROR] --sample is required" >&2
    exit 1
fi

# Check seqtk
if ! command -v seqtk &>/dev/null; then
    echo "[ERROR] seqtk not found. Install: conda install -c bioconda seqtk" >&2
    exit 1
fi

FASTQ_DIR="${DS_RAW}/fastq/${SAMPLE}"
SUBSET_FASTQ="${DS_RAW}/fastq_subset/${SAMPLE}_subset"
SUBSET_OUT="${DS_RAW}/cellranger_out/subset_test"

mkdir -p "${SUBSET_FASTQ}" "${SUBSET_OUT}"

echo ""
echo "=========================================="
echo "Subset Test: ${SAMPLE}"
echo "N reads:     ${N_READS}"
echo "=========================================="

# Subsample each FASTQ pair
for R1 in "${FASTQ_DIR}"/*_R1_001.fastq.gz; do
    R2="${R1/_R1_/_R2_}"
    BASE=$(basename "${R1}")
    OUT_R1="${SUBSET_FASTQ}/${BASE}"
    OUT_R2="${SUBSET_FASTQ}/$(basename "${R2}")"

    if [ -f "${OUT_R1}" ] && [ -f "${OUT_R2}" ]; then
        echo "[info] Subset already exists: ${BASE}, skipping"
        continue
    fi

    echo "[info] Subsampling: ${BASE} → ${N_READS} reads"
    seqtk sample -s${SEED} "${R1}" ${N_READS} | gzip > "${OUT_R1}"
    seqtk sample -s${SEED} "${R2}" ${N_READS} | gzip > "${OUT_R2}"
done

echo "[info] Subset FASTQ: ${SUBSET_FASTQ}"
ls -lh "${SUBSET_FASTQ}"

# Run Cell Ranger on subset
echo ""
echo "[info] Running Cell Ranger count on subset..."

cd "${SUBSET_OUT}"
rm -rf "${SAMPLE}_subset"

PREFIXES=$(get_sample_prefixes "${SUBSET_FASTQ}")

"${CELLRANGER_PATH}/cellranger" count \
    --id="${SAMPLE}_subset" \
    --fastqs="${SUBSET_FASTQ}" \
    --sample="${PREFIXES}" \
    --transcriptome="${CELLRANGER_REF}" \
    --create-bam=false \
    --chemistry=auto \
    --include-introns=true \
    --localcores="${CORES}" \
    --localmem="${MEM_GB}"

EXIT_CODE=$?

echo ""
if [ ${EXIT_CODE} -eq 0 ]; then
    echo "✅ Subset test PASSED"
    echo ""
    echo "--- Quick metrics ---"
    extract_qc_metrics "${SUBSET_OUT}/${SAMPLE}_subset/outs/metrics_summary.csv" "${SAMPLE}_subset"
else
    echo "❌ Subset test FAILED (exit code: ${EXIT_CODE})"
fi

exit ${EXIT_CODE}
