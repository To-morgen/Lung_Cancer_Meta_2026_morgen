# ============================================================================
# gates.smk — Manual review checkpoints
#
# Usage:
#   1. Run pipeline up to a gate:  snakemake phase05_done --cores 16
#   2. Review results in Cursor
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

rule gate_annotation_review:
    """MANUAL GATE: Review annotation before DE.
    To approve:  touch {output.gate}"""
    input:
        os.path.join(RESULTS, "06_annotate", ".phase06_done")
    output:
        gate = os.path.join(RESULTS, "06_annotate", ".gate_approved")
    shell:
        """
        echo ""
        echo "╔═══════════════════════════════════════════════════════════╗"
        echo "║  🔍 MANUAL GATE: Annotation Review Required              ║"
        echo "║                                                           ║"
        echo "║  Review:                                                  ║"
        echo "║    {RESULTS}/06_annotate/plots/                           ║"
        echo "║    configs/annotation/celltype_mapping.csv                ║"
        echo "║                                                           ║"
        echo "║  When approved, run:                                      ║"
        echo "║    touch {output.gate}                                    ║"
        echo "║                                                           ║"
        echo "║  Then re-run snakemake to continue.                       ║"
        echo "╚═══════════════════════════════════════════════════════════╝"
        echo ""
        exit 1
        """
