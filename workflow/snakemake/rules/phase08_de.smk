rule phase08_de:
    """Phase 08: Pseudobulk DE → Enrichment → Visualization"""
    input:
        os.path.join(RESULTS, "06_annotate", ".gate_approved")
    output:
        sentinel = os.path.join(RESULTS, "08_de", ".phase08_done")
    params:
        project_root = PROJECT_ROOT,
        env_vars     = ENV_VARS,
        rscript      = RSCRIPT
    log:
        os.path.join(LOGS, "phase08_de.log")
    shell:
        """
        cd {params.project_root}
        {params.env_vars} {params.rscript} workflow/scrna/08_de/run_de_pipeline.R \
            2>&1 | tee {log}
        touch {output.sentinel}
        """
