#!/bin/bash
# ============================================================================
# 02_cellranger_qc_summary.sh — Generate QC summary from Cell Ranger outputs
#
# Usage:
#   bash workflow/scrna/01_alignment/02_cellranger_qc_summary.sh \
#       --output-subdir full \
#       --samples "A1_1 A1_2 FL_1 FL_2 mc_1 mc_2"
#
# Or auto-detect all completed samples:
#   bash workflow/scrna/01_alignment/02_cellranger_qc_summary.sh --output-subdir full --auto
#
# Outputs:
#   - Console: formatted QC table
#   - File:    {DS_RAW}/cellranger_out/{subdir}/qc_summary.csv
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../functions/cellranger_utils.sh"

OUTPUT_SUBDIR="full"
SAMPLES=""
AUTO_DETECT=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --output-subdir) OUTPUT_SUBDIR="$2"; shift 2 ;;
        --samples)       SAMPLES="$2";       shift 2 ;;
        --auto)          AUTO_DETECT=true;   shift ;;
        -h|--help)       head -15 "$0" | grep '^#' | sed 's/^# *//'; exit 0 ;;
        *)               echo "[ERROR] Unknown: $1" >&2; exit 1 ;;
    esac
done

load_project_env
OUTPUT_DIR="${DS_RAW}/cellranger_out/${OUTPUT_SUBDIR}"

if [ "${AUTO_DETECT}" = true ]; then
    SAMPLES=$(ls -d "${OUTPUT_DIR}"/*/outs/metrics_summary.csv 2>/dev/null \
        | xargs -I{} dirname {} | xargs -I{} dirname {} | xargs -I{} basename {} \
        | tr '\n' ' ')
fi

if [ -z "${SAMPLES}" ]; then
    echo "[ERROR] No samples found/specified" >&2
    exit 1
fi

SUMMARY_CSV="${OUTPUT_DIR}/qc_summary.csv"
HEADER="Sample,Estimated Number of Cells,Mean Reads per Cell,Median Genes per Cell,Median UMI Counts per Cell,Sequencing Saturation,Reads Mapped Confidently to Transcriptome,Fraction Reads in Cells,Total Genes Detected,Valid Barcodes"

echo "${HEADER}" > "${SUMMARY_CSV}"

echo ""
echo "========== Cell Ranger QC Summary =========="
echo "Source: ${OUTPUT_DIR}"
echo ""
printf "%-8s %10s %10s %12s %12s %10s %10s %10s\n" \
    "Sample" "Cells" "MeanReads" "MedGenes" "MedUMI" "Satur." "MapTrans" "FracCell"
echo "--------------------------------------------------------------------------------------------"

for SAMPLE in ${SAMPLES}; do
    CSV="${OUTPUT_DIR}/${SAMPLE}/outs/metrics_summary.csv"
    ROW=$(extract_qc_metrics "${CSV}" "${SAMPLE}")
    echo "${ROW}" >> "${SUMMARY_CSV}"

    # Pretty print
    python3 -c "
row = '${ROW}'.split(',')
if len(row) >= 10:
    print(f'{row[0]:<8s} {row[1]:>10s} {row[2]:>10s} {row[3]:>12s} {row[4]:>12s} {row[5]:>10s} {row[6]:>10s} {row[7]:>10s}')
else:
    print(f'{row[0]:<8s} INCOMPLETE')
" 2>/dev/null || echo "${SAMPLE}  (parse error)"
done

echo ""
echo "Summary saved: ${SUMMARY_CSV}"
echo ""

# === Flag checks ===
echo "========== Auto-QC Flags =========="
python3 << PYEOF
import csv
with open('${SUMMARY_CSV}') as f:
    reader = csv.DictReader(f)
    flags = []
    for row in reader:
        s = row['Sample']
        try:
            frac = float(row.get('Fraction Reads in Cells','0').replace('%','').replace('"',''))
            if frac < 80:
                flags.append(f"  ⚠️  {s}: Fraction Reads in Cells = {frac}% (<80%)")
        except: pass
        try:
            vb = float(row.get('Valid Barcodes','0').replace('%','').replace('"',''))
            if vb < 90:
                flags.append(f"  ⚠️  {s}: Valid Barcodes = {vb}% (<90%)")
        except: pass
        try:
            mg = int(row.get('Median Genes per Cell','0').replace(',','').replace('"',''))
            if mg < 1000:
                flags.append(f"  ⚠️  {s}: Median Genes = {mg} (<1000)")
        except: pass
    if flags:
        print('\n'.join(flags))
    else:
        print("  ✅ All samples passed auto-QC checks")
PYEOF
echo ""
