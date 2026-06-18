rule phase08_de:
    """Phase 08: Pseudobulk DE → Enrichment → Visualization"""
    input:
        gate = os.path.join(RESULTS, "06_annotate", ".gate_approved"),
        final_object = os.path.join(RESULTS, "06_annotate", "objects", "seurat_annotated_final.rds")
    output:
        sentinel = os.path.join(RESULTS, "08_de", ".phase08_done")
    params:
        project_root = PROJECT_ROOT,
        env_vars     = ENV_VARS,
        de_config    = DE_CONTRASTS_CONFIG,
        rscript      = RSCRIPT
    log:
        os.path.join(LOGS, "phase08_de.log")
    shell:
        """
        cd {params.project_root}
        {params.env_vars} DE_CONTRASTS_CONFIG="{params.de_config}" {params.rscript} workflow/scrna/08_de/run_de_pipeline.R \
            2>&1 | tee {log}
        touch {output.sentinel}
        """
