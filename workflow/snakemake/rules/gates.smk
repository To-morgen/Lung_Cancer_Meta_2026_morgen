# ============================================================================
# gates.smk — Manual review checkpoints
#
# Usage:
#   1. Run pipeline up to a gate:  snakemake phase05_done --cores 16
#   2. Review results in your editor/browser
#   3. Approve:  touch results/scrna/{ds_prefix}/05_cluster/.gate_approved
#   4. Continue:  snakemake phase08_done --cores 16
# ============================================================================

rule gate_cluster_review:
    """MANUAL GATE: Review clustering before annotation.
    To approve:  touch {output.gate}"""
    input:
        os.path.join(RESULTS, "05_cluster", ".phase05_done")
    output:
        gate = os.path.join(RESULTS, "05_cluster", ".gate_approved")
    shell:
        """
        echo ""
        echo "╔═══════════════════════════════════════════════════════════╗"
        echo "║  🔍 MANUAL GATE: Cluster Review Required                 ║"
        echo "║                                                           ║"
        echo "║  Review plots in:                                         ║"
        echo "║    {RESULTS}/05_cluster/plots/                            ║"
        echo "║                                                           ║"
        echo "║  When approved, run:                                      ║"
        echo "║    touch {output.gate}                                    ║"
        echo "║                                                           ║"
        echo "║  Then re-run snakemake to continue.                       ║"
        echo "╚═══════════════════════════════════════════════════════════╝"
        echo ""
        exit 1
        """

rule gate_annotation_mapping_review:
    """MANUAL GATE: Review auto annotation and approve cell-type mapping."""
    input:
        os.path.join(RESULTS, "06_annotate", ".phase06_auto_done")
    output:
        gate = os.path.join(RESULTS, "06_annotate", ".manual_annotation_approved")
    shell:
        """
        echo ""
        echo "╔═══════════════════════════════════════════════════════════╗"
        echo "║  🔍 MANUAL GATE: Annotation Mapping Review Required      ║"
        echo "║                                                           ║"
        echo "║  Review:                                                  ║"
        echo "║    {RESULTS}/06_annotate/plots/                           ║"
        echo "║    {RESULTS}/06_annotate/markers/                         ║"
        echo "║    {RESULTS}/06_annotate/reports/singler_cluster_annotation.csv"
        echo "║                                                           ║"
        echo "║  Fill or confirm mapping:                                 ║"
        echo "║    configs/annotation/celltype_mapping_{DS_PREFIX}.csv    ║"
        echo "║    or fallback configs/annotation/celltype_mapping.csv    ║"
        echo "║                                                           ║"
        echo "║  When mapping is approved, run:                           ║"
        echo "║    touch {output.gate}                                    ║"
        echo "║                                                           ║"
        echo "║  Then re-run snakemake to apply manual annotation.        ║"
        echo "╚═══════════════════════════════════════════════════════════╝"
        echo ""
        exit 1
        """

rule gate_annotation_review:
    """MANUAL GATE: Review final annotation before DE.
    To approve:  touch {output.gate}"""
    input:
        sentinel = os.path.join(RESULTS, "06_annotate", ".phase06_manual_done"),
        final_object = os.path.join(RESULTS, "06_annotate", "objects", "seurat_annotated_final.rds")
    output:
        gate = os.path.join(RESULTS, "06_annotate", ".gate_approved")
    shell:
        """
        echo ""
        echo "╔═══════════════════════════════════════════════════════════╗"
        echo "║  🔍 MANUAL GATE: Final Annotation Review Required        ║"
        echo "║                                                           ║"
        echo "║  Review final annotation outputs:                         ║"
        echo "║    {RESULTS}/06_annotate/plots/                           ║"
        echo "║    {RESULTS}/06_annotate/reports/                         ║"
        echo "║    {RESULTS}/06_annotate/objects/seurat_annotated_final.rds"
        echo "║                                                           ║"
        echo "║  When approved, run:                                      ║"
        echo "║    touch {output.gate}                                    ║"
        echo "║                                                           ║"
        echo "║  Then re-run snakemake to continue.                       ║"
        echo "╚═══════════════════════════════════════════════════════════╝"
        echo ""
        exit 1
        """
