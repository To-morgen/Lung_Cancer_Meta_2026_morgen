rule phase06_annotate:
    """Phase 06: FindAllMarkers → SingleR → Manual Annotation"""
    input:
        os.path.join(RESULTS, "05_cluster", ".gate_approved")
    output:
        sentinel = os.path.join(RESULTS, "06_annotate", ".phase06_done")
    params:
        project_root = PROJECT_ROOT,
        env_vars     = ENV_VARS,
        rscript      = RSCRIPT
    log:
        os.path.join(LOGS, "phase06_annotate.log")
    shell:
        """
        cd {params.project_root}
        {params.env_vars} {params.rscript} workflow/scrna/06_annotate/run_annotation_pipeline.R \
            2>&1 | tee {log}
        touch {output.sentinel}
        """
