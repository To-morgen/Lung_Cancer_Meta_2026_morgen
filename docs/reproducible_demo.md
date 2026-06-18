# Reproducible Demo

Milestone 2 adds a lightweight smoke demo that can run from a clean public clone
without private data.

The demo does **not** perform biological analysis. It verifies that the public
repository can:

1. generate a tiny 10x-like input structure,
2. load a Snakemake workflow,
3. validate the expected dataset input contract, and
4. write a reproducible local sentinel/report under ignored output paths.

## Run the demo

```bash
Rscript examples/toy_lung_2grp/generate_toy_inputs.R
snakemake -s workflow/datasets/smoke_lung_2grp/Snakefile -n --cores 1
snakemake -s workflow/datasets/smoke_lung_2grp/Snakefile --cores 1
```

Or use the helper:

```bash
bash scripts/verify_public_demo.sh

# If Snakemake is not on PATH:
SNAKEMAKE_BIN=/path/to/snakemake bash scripts/verify_public_demo.sh
```

## Expected outputs

```text
data/demo/toy_lung_2grp/                         # generated inputs, ignored
results/demo/smoke_lung_2grp/input_contract_report.tsv
results/demo/smoke_lung_2grp/.smoke_demo_done
logs/demo/smoke_lung_2grp/smoke_validate_inputs.log
```

These outputs are ignored by Git and can be safely regenerated.

## What this proves

This smoke demo proves public reproducibility at the repository-contract level:
configuration paths, generated input layout, Snakemake target resolution, and
basic execution all work without private datasets.

Full biological phases such as Seurat QC, CNV inference, DESeq2, and CellChat
remain real-data workflows and require appropriate local data and environments.
