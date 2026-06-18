rule phase06_annotate_auto:
    """Phase 06a: FindAllMarkers → SingleR auto annotation"""
    input:
        os.path.join(RESULTS, "05_cluster", ".gate_approved")
    output:
        sentinel = os.path.join(RESULTS, "06_annotate", ".phase06_auto_done")
    params:
        project_root = PROJECT_ROOT,
        env_vars     = ENV_VARS,
        rscript      = RSCRIPT
    log:
        os.path.join(LOGS, "phase06_annotate_auto.log")
    shell:
        """
        cd {params.project_root}
        {params.env_vars} {params.rscript} workflow/scrna/06_annotate/run_annotation_pipeline.R \
            2>&1 | tee {log}
        touch {output.sentinel}
        """

rule phase06_annotate_manual:
    """Phase 06b: Apply reviewed cell-type mapping"""
    input:
        gate = os.path.join(RESULTS, "06_annotate", ".manual_annotation_approved")
    output:
        sentinel = os.path.join(RESULTS, "06_annotate", ".phase06_manual_done"),
        legacy_sentinel = os.path.join(RESULTS, "06_annotate", ".phase06_done"),
        final_object = os.path.join(RESULTS, "06_annotate", "objects", "seurat_annotated_final.rds")
    params:
        project_root = PROJECT_ROOT,
        env_vars     = ENV_VARS,
        rscript      = RSCRIPT
    log:
        os.path.join(LOGS, "phase06_annotate_manual.log")
    shell:
        """
        cd {params.project_root}
        {params.env_vars} {params.rscript} workflow/scrna/06_annotate/run_annotation_pipeline.R --manual \
            2>&1 | tee {log}
        touch {output.sentinel} {output.legacy_sentinel}
        """
