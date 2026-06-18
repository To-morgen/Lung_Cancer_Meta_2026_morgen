rule phase02_qc:
    """Phase 02: SoupX → Create Seurat → MAD QC → Doublet Removal → Viz"""
    input:
        [os.path.join(RESULTS, "01_alignment", ".phase01_done")]
        if config.get("run_phases", {}).get("phase01_alignment", False)
        else (
            expand(
                os.path.join(config["cellranger_out"], "{sample}", "outs",
                             "filtered_feature_bc_matrix"),
                sample=SAMPLES
            ) if "cellranger_out" in config and config["cellranger_out"] else []
        )
    output:
        sentinel = os.path.join(RESULTS, "02_qc", ".phase02_done")
    params:
        project_root = PROJECT_ROOT,
        env_vars     = ENV_VARS,
        rscript      = RSCRIPT
    log:
        os.path.join(LOGS, "phase02_qc.log")
    shell:
        """
        cd {params.project_root}
        {params.env_vars} {params.rscript} workflow/scrna/02_qc/run_qc_pipeline.R \
            2>&1 | tee {log}
        touch {output.sentinel}
        """
