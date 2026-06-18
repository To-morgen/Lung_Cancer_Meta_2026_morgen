rule phase09_cellchat_prepare:
    """Phase 09a: Prepare per-group CellChat objects."""
    input:
        gate = os.path.join(RESULTS, "06_annotate", ".gate_approved"),
        phase06_done = os.path.join(RESULTS, "06_annotate", ".phase06_manual_done"),
        final_object = os.path.join(RESULTS, "06_annotate", "objects", "seurat_annotated_final.rds"),
        de_config = DE_CONTRASTS_CONFIG
    output:
        sentinel = os.path.join(RESULTS, "09_cellchat", ".phase09_cellchat_prepare_done"),
        design = os.path.join(RESULTS, "09_cellchat", "objects", "design.rds"),
        cellchat_list = os.path.join(RESULTS, "09_cellchat", "objects", "cellchat_list.rds")
    threads: 4
    resources:
        mem_mb = 64000
    params:
        project_root = PROJECT_ROOT,
        cellchat_dir = os.path.join(PROJECT_ROOT, "modules", "cellchat"),
        rscript = RSCRIPT
    log:
        os.path.join(LOGS, "phase09_cellchat_prepare.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {log})
        cd {params.cellchat_dir}
        LUNGMETA_ROOT="{params.project_root}" DATASET_ID="{DS_PREFIX}" DE_CONTRASTS_CONFIG="{input.de_config}" \
          {params.rscript} scripts/01_prepare_cellchat.R \
          2>&1 | tee {log}
        touch {output.sentinel}
        """


rule phase09_cellchat_run:
    """Phase 09b: Run CellChat inference per group."""
    input:
        prepare_done = os.path.join(RESULTS, "09_cellchat", ".phase09_cellchat_prepare_done"),
        design = os.path.join(RESULTS, "09_cellchat", "objects", "design.rds"),
        cellchat_list = os.path.join(RESULTS, "09_cellchat", "objects", "cellchat_list.rds")
    output:
        sentinel = os.path.join(RESULTS, "09_cellchat", ".phase09_cellchat_run_done"),
        inferred = os.path.join(RESULTS, "09_cellchat", "objects", "cellchat_list_inferred.rds")
    threads: 16
    resources:
        mem_mb = 128000
    params:
        project_root = PROJECT_ROOT,
        cellchat_dir = os.path.join(PROJECT_ROOT, "modules", "cellchat"),
        rscript = RSCRIPT
    log:
        os.path.join(LOGS, "phase09_cellchat_run.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {log})
        cd {params.cellchat_dir}
        LUNGMETA_ROOT="{params.project_root}" DATASET_ID="{DS_PREFIX}" \
          {params.rscript} scripts/02_run_cellchat.R \
          2>&1 | tee {log}
        touch {output.sentinel}
        """


rule phase09_cellchat_compare:
    """Phase 09c: Compare CellChat groups."""
    input:
        run_done = os.path.join(RESULTS, "09_cellchat", ".phase09_cellchat_run_done"),
        inferred = os.path.join(RESULTS, "09_cellchat", "objects", "cellchat_list_inferred.rds"),
        design = os.path.join(RESULTS, "09_cellchat", "objects", "design.rds")
    output:
        sentinel = os.path.join(RESULTS, "09_cellchat", ".phase09_cellchat_compare_done"),
        merged = os.path.join(RESULTS, "09_cellchat", "objects", "cellchat_merged.rds")
    threads: 4
    resources:
        mem_mb = 128000
    params:
        project_root = PROJECT_ROOT,
        cellchat_dir = os.path.join(PROJECT_ROOT, "modules", "cellchat"),
        rscript = RSCRIPT
    log:
        os.path.join(LOGS, "phase09_cellchat_compare.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {log})
        cd {params.cellchat_dir}
        LUNGMETA_ROOT="{params.project_root}" DATASET_ID="{DS_PREFIX}" \
          {params.rscript} scripts/03_compare_groups.R \
          2>&1 | tee {log}
        touch {output.sentinel}
        """


rule phase09_cellchat_deep_dive:
    """Phase 09d: Run CellChat pathway and pair deep-dive plots."""
    input:
        compare_done = os.path.join(RESULTS, "09_cellchat", ".phase09_cellchat_compare_done"),
        inferred = os.path.join(RESULTS, "09_cellchat", "objects", "cellchat_list_inferred.rds"),
        merged = os.path.join(RESULTS, "09_cellchat", "objects", "cellchat_merged.rds"),
        design = os.path.join(RESULTS, "09_cellchat", "objects", "design.rds")
    output:
        sentinel = os.path.join(RESULTS, "09_cellchat", ".phase09_cellchat_deep_dive_done")
    threads: 4
    resources:
        mem_mb = 128000
    params:
        project_root = PROJECT_ROOT,
        cellchat_dir = os.path.join(PROJECT_ROOT, "modules", "cellchat"),
        rscript = RSCRIPT
    log:
        os.path.join(LOGS, "phase09_cellchat_deep_dive.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {log})
        cd {params.cellchat_dir}
        LUNGMETA_ROOT="{params.project_root}" DATASET_ID="{DS_PREFIX}" \
          {params.rscript} scripts/04_deep_dive.R \
          2>&1 | tee {log}
        touch {output.sentinel}
        """


rule phase09_cellchat_summary:
    """Phase 09e: Summarize generated CellChat CSV reports."""
    input:
        run_done = os.path.join(RESULTS, "09_cellchat", ".phase09_cellchat_run_done"),
        compare_done = os.path.join(RESULTS, "09_cellchat", ".phase09_cellchat_compare_done"),
        deep_dive_done = os.path.join(RESULTS, "09_cellchat", ".phase09_cellchat_deep_dive_done"),
        inference_summary = os.path.join(RESULTS, "09_cellchat", "qc", "inference_summary.csv"),
        pathway_presence = os.path.join(RESULTS, "09_cellchat", "comparison", "reports", "pathway_presence_matrix.csv")
    output:
        sentinel = os.path.join(RESULTS, "09_cellchat", ".phase09_cellchat_summary_done"),
        summary_index = os.path.join(RESULTS, "09_cellchat", "summary", "reports", "00_summary_index.csv")
    threads: 1
    resources:
        mem_mb = 16000
    params:
        project_root = PROJECT_ROOT,
        cellchat_dir = os.path.join(PROJECT_ROOT, "modules", "cellchat"),
        rscript = RSCRIPT
    log:
        os.path.join(LOGS, "phase09_cellchat_summary.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {log})
        cd {params.cellchat_dir}
        LUNGMETA_ROOT="{params.project_root}" DATASET_ID="{DS_PREFIX}" \
          {params.rscript} scripts/05_summarize_reports.R \
          2>&1 | tee {log}
        touch {output.sentinel}
        """
