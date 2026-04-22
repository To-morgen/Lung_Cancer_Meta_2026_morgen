rule phase04_integrate:
    """Phase 04: Harmony batch correction"""
    input:
        os.path.join(RESULTS, "03_normalize", ".phase03_done")
    output:
        sentinel = os.path.join(RESULTS, "04_integrate", ".phase04_done")
    params:
        project_root = PROJECT_ROOT,
        env_vars     = ENV_VARS,
        rscript      = RSCRIPT
    log:
        os.path.join(LOGS, "phase04_integrate.log")
    shell:
        """
        cd {params.project_root}
        {params.env_vars} {params.rscript} workflow/scrna/04_integrate/run_integration_pipeline.R \
            2>&1 | tee {log}
        touch {output.sentinel}
        """
