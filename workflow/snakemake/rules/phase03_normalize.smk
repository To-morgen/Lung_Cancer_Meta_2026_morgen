rule phase03_normalize:
    """Phase 03: Merge → Cell Cycle Score → SCTransform → PCA → CC Assessment"""
    input:
        os.path.join(RESULTS, "02_qc", ".phase02_done")
    output:
        sentinel = os.path.join(RESULTS, "03_normalize", ".phase03_done")
    params:
        project_root = PROJECT_ROOT,
        env_vars     = ENV_VARS,
        rscript      = RSCRIPT
    log:
        os.path.join(LOGS, "phase03_normalize.log")
    shell:
        """
        cd {params.project_root}
        {params.env_vars} {params.rscript} workflow/scrna/03_normalize/run_normalize_pipeline.R \
            2>&1 | tee {log}
        touch {output.sentinel}
        """
