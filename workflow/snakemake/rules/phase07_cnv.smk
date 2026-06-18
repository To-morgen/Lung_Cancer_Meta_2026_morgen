rule phase07_scevan:
    """Phase 07a: Run full SCEVAN CNV screening in isolated CNV module."""
    input:
        final_object = os.path.join(RESULTS, "06_annotate", "objects", "seurat_annotated_final.rds"),
        cnv_config = os.path.join(PROJECT_ROOT, "configs", "annotation", f"cnv_targets_{DS_PREFIX}.yaml")
    output:
        sentinel = os.path.join(RESULTS, "07_cnv", ".phase07_scevan_done"),
        scevan_object = os.path.join(RESULTS, "07_cnv", "scevan", "seurat_with_scevan.rds"),
        cluster_report = os.path.join(RESULTS, "07_cnv", "reports", "scevan_tumor_fraction_by_cluster.csv")
    threads: 8
    params:
        project_root = PROJECT_ROOT,
        cnv_dir = os.path.join(PROJECT_ROOT, "modules", "cnv"),
        rscript = RSCRIPT
    log:
        os.path.join(LOGS, "phase07_scevan.log")
    shell:
        """
        mkdir -p $(dirname {log})
        cd {params.cnv_dir}
        LUNGMETA_ROOT="{params.project_root}" CNV_CONFIG="{input.cnv_config}" SCEVAN_CORES="{threads}" \
          {params.rscript} scripts/01_scevan.R "{input.cnv_config}" \
          2>&1 | tee {log}
        LUNGMETA_ROOT="{params.project_root}" CNV_CONFIG="{input.cnv_config}" SCEVAN_CORES="{threads}" \
          {params.rscript} scripts/02_validate_scevan.R "{input.cnv_config}" \
          2>&1 | tee -a {log}
        touch {output.sentinel}
        """


rule phase07_infercnv_prep:
    """Phase 07b: Prepare dataset-specific inferCNV inputs."""
    input:
        final_object = os.path.join(RESULTS, "06_annotate", "objects", "seurat_annotated_final.rds"),
        cnv_config = os.path.join(PROJECT_ROOT, "configs", "annotation", f"cnv_targets_{DS_PREFIX}.yaml")
    output:
        sentinel = os.path.join(RESULTS, "07_cnv", ".phase07_infercnv_prep_done"),
        counts = os.path.join(RESULTS, "07_cnv", "infercnv", "inputs", "counts_sparse.rds"),
        annotations = os.path.join(RESULTS, "07_cnv", "infercnv", "inputs", "cell_annotations.tsv"),
        reference_groups = os.path.join(RESULTS, "07_cnv", "infercnv", "inputs", "reference_groups.txt"),
        dataset_info = os.path.join(RESULTS, "07_cnv", "infercnv", "inputs", "dataset_info.yaml"),
        prep_summary = os.path.join(RESULTS, "07_cnv", "infercnv", "inputs", "prep_summary.csv")
    params:
        project_root = PROJECT_ROOT,
        cnv_dir = os.path.join(PROJECT_ROOT, "modules", "cnv"),
        rscript = RSCRIPT
    log:
        os.path.join(LOGS, "phase07_infercnv_prep.log")
    shell:
        """
        mkdir -p $(dirname {log})
        cd {params.cnv_dir}
        LUNGMETA_ROOT="{params.project_root}" CNV_CONFIG="{input.cnv_config}" \
          {params.rscript} scripts/06_prep_infercnv_inputs.R "{input.cnv_config}" \
          2>&1 | tee {log}
        touch {output.sentinel}
        """


rule phase07_infercnv:
    """Phase 07c: Run inferCNV from prepared dataset-specific inputs."""
    input:
        prep_done = os.path.join(RESULTS, "07_cnv", ".phase07_infercnv_prep_done"),
        cnv_config = os.path.join(PROJECT_ROOT, "configs", "annotation", f"cnv_targets_{DS_PREFIX}.yaml"),
        counts = os.path.join(RESULTS, "07_cnv", "infercnv", "inputs", "counts_sparse.rds"),
        annotations = os.path.join(RESULTS, "07_cnv", "infercnv", "inputs", "cell_annotations.tsv"),
        reference_groups = os.path.join(RESULTS, "07_cnv", "infercnv", "inputs", "reference_groups.txt"),
        dataset_info = os.path.join(RESULTS, "07_cnv", "infercnv", "inputs", "dataset_info.yaml")
    output:
        sentinel = os.path.join(RESULTS, "07_cnv", ".phase07_infercnv_done"),
        infercnv_object = os.path.join(RESULTS, "07_cnv", "infercnv", "run", "infercnv_obj_final.rds")
    threads: 16
    params:
        project_root = PROJECT_ROOT,
        cnv_dir = os.path.join(PROJECT_ROOT, "modules", "cnv"),
        rscript = RSCRIPT
    log:
        os.path.join(LOGS, "phase07_infercnv.log")
    shell:
        """
        mkdir -p $(dirname {log})
        cd {params.cnv_dir}
        LUNGMETA_ROOT="{params.project_root}" CNV_CONFIG="{input.cnv_config}" INFERCNV_THREADS="{threads}" \
          {params.rscript} scripts/02_infercnv.R "{input.cnv_config}" \
          2>&1 | tee {log}
        touch {output.sentinel}
        """


rule phase07_infercnv_score:
    """Phase 07d: Compute per-cell inferCNV burden scores and diagnostic plots."""
    input:
        infercnv_done = os.path.join(RESULTS, "07_cnv", ".phase07_infercnv_done"),
        final_object = os.path.join(RESULTS, "06_annotate", "objects", "seurat_annotated_final.rds"),
        cnv_config = os.path.join(PROJECT_ROOT, "configs", "annotation", f"cnv_targets_{DS_PREFIX}.yaml"),
        infercnv_object = os.path.join(RESULTS, "07_cnv", "infercnv", "run", "infercnv_obj_final.rds"),
        annotations = os.path.join(RESULTS, "07_cnv", "infercnv", "inputs", "cell_annotations.tsv")
    output:
        sentinel = os.path.join(RESULTS, "07_cnv", ".phase07_infercnv_score_done"),
        per_cell_scores = os.path.join(RESULTS, "07_cnv", "infercnv", "scoring", "reports", "infercnv_per_cell_scores.csv"),
        by_cluster = os.path.join(RESULTS, "07_cnv", "infercnv", "scoring", "reports", "infercnv_score_by_cluster.csv"),
        overview_plot = os.path.join(RESULTS, "07_cnv", "infercnv", "scoring", "plots", "01_infercnv_score_overview.pdf")
    params:
        project_root = PROJECT_ROOT,
        cnv_dir = os.path.join(PROJECT_ROOT, "modules", "cnv"),
        rscript = RSCRIPT
    log:
        os.path.join(LOGS, "phase07_infercnv_score.log")
    shell:
        """
        mkdir -p $(dirname {log})
        cd {params.cnv_dir}
        LUNGMETA_ROOT="{params.project_root}" CNV_CONFIG="{input.cnv_config}" \
          {params.rscript} scripts/07_score_infercnv.R "{input.cnv_config}" \
          2>&1 | tee {log}
        touch {output.sentinel}
        """
