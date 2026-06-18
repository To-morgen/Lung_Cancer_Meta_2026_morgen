# Toy Lung 2-Group Demo

This demo creates a tiny 10x-like directory structure that is sufficient for
public smoke tests. It is not intended for biological interpretation.

## Generate inputs

```bash
Rscript examples/toy_lung_2grp/generate_toy_inputs.R
```

This writes small Matrix Market files under:

```text
data/demo/toy_lung_2grp/cellranger_out/{sample}/outs/filtered_feature_bc_matrix/
```

The generated files are ignored by Git.

## Run smoke workflow

```bash
snakemake -s workflow/datasets/smoke_lung_2grp/Snakefile -n --cores 1
snakemake -s workflow/datasets/smoke_lung_2grp/Snakefile --cores 1
```

The smoke workflow checks that the expected input contract exists and writes a
small report under `results/demo/smoke_lung_2grp/`.
