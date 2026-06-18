rule phase01_alignment:
    """Phase 01: Cell Ranger count → QC summary.
    Produces cellranger_out/{output_subdir}/{sample}/outs/ which phase02 consumes.
    """
    input:
        # FASTQ are not modelled as Snakemake file inputs (large, private).
        # Cellranger expects DS_RAW/fastq/{sample}/*.fastq.gz at runtime.
    output:
        sentinel = os.path.join(RESULTS, "01_alignment", ".phase01_done")
    params:
        project_root = PROJECT_ROOT,
        dataset_id   = config["dataset_id"],
        samples      = " ".join(SAMPLES),
        output_subdir = config.get("phase01", {}).get("output_subdir", "full"),
        cores   = config.get("phase01", {}).get("cores", 32),
        mem_gb  = config.get("phase01", {}).get("mem_gb", 128),
        bam     = str(config.get("phase01", {}).get("create_bam", True)).lower()
    log:
        os.path.join(LOGS, "phase01_alignment.log")
    shell:
        """
        mkdir -p {RESULTS}/01_alignment $(dirname {log})

        # 1) Cell Ranger count (skips already-completed samples)
        bash workflow/scrna/01_alignment/01_cellranger_count.sh \
            --dataset {params.dataset_id} \
            --samples "{params.samples}" \
            --cores {params.cores} \
            --mem {params.mem_gb} \
            --bam {params.bam} \
            --output-subdir {params.output_subdir} \
            2>&1 | tee {log}

        # 2) QC summary
        bash workflow/scrna/01_alignment/02_cellranger_qc_summary.sh \
            --output-subdir {params.output_subdir} \
            --samples "{params.samples}" \
            2>&1 | tee -a {log}

        touch {output.sentinel}
        """
