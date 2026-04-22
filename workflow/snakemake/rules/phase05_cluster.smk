rule phase05_cluster:
    """Phase 05: UMAP + Leiden Clustering + Cluster QC"""
    input:
        os.path.join(RESULTS, "04_integrate", ".phase04_done")
    output:
        sentinel = os.path.join(RESULTS, "05_cluster", ".phase05_done")
    params:
        project_root = PROJECT_ROOT,
        env_vars     = ENV_VARS,
        rscript      = RSCRIPT
    log:
        os.path.join(LOGS, "phase05_cluster.log")
    shell:
        """
        cd {params.project_root}
        {params.env_vars} {params.rscript} workflow/scrna/05_cluster/run_cluster_pipeline.R \
            2>&1 | tee {log}
        touch {output.sentinel}
        """
