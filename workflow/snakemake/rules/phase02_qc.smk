rule phase02_qc:
    """Phase 02: SoupX → Create Seurat → MAD QC → Doublet → Viz"""
    output:
        sentinel = os.path.join(RESULTS, "02_qc", ".phase02_done")
    params:
        project_root = PROJECT_ROOT,
        env_vars     = ENV_VARS
    log:
        os.path.join(LOGS, "phase02_qc.log")
    shell:
        """
        cd {params.project_root}
        {params.env_vars} Rscript workflow/scrna/02_qc/run_qc_pipeline.R \
            2>&1 | tee {log}
        touch {output.sentinel}
        """
