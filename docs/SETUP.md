# Setup and Governance

This project separates reusable pipeline code from private data and generated
analysis outputs.

## Storage policy

| Category | Track in Git | Keep out of Git |
|---|---|---|
| Pipeline logic | `workflow/`, `modules/`, `scripts/` | ad hoc local scratch code |
| Configuration | sanitized YAML/CSV templates | private paths, sample manifests, unpublished decisions |
| Environment | `renv.lock`, module `renv.lock`, activation files | installed libraries, conda env directories |
| Data | small schemas/templates only | FASTQ, BAM, MTX, RDS, H5, H5AD, CellChat/inferCNV objects |
| Outputs | none by default | `results/`, `figures/`, `logs/`, Snakemake runtime state |
| Notes | public docs only | handoff files, working memory, local agent/tooling state |

Private dataset config files should live under `configs_private/`, which is
ignored by Git. Public configs must use relative paths, placeholders, or example
values only.

## Environment

The main project and each heavy downstream module use separate R environments.

```bash
# Main project
Rscript -e 'renv::restore()'

# CNV module
(cd modules/cnv && Rscript -e 'renv::restore()')

# CellChat module
(cd modules/cellchat && Rscript -e 'renv::restore()')
```

Snakemake should be installed in your preferred workflow environment. A typical
local execution pattern is:

```bash
snakemake -s workflow/datasets/<dataset>/Snakefile -n --cores 1
snakemake -s workflow/datasets/<dataset>/Snakefile --cores 8
```

## Local environment file

Copy the example and edit paths for your machine:

```bash
cp .env.sh.example .env.sh
```

`.env.sh` is ignored by Git. Do not commit machine-specific paths, credentials,
or unpublished dataset locations.

## Dataset entry points

A dataset workflow contains:

```text
workflow/datasets/<dataset>/
  Snakefile
  config.yaml
```

The Snakemake `config.yaml` should contain lightweight orchestration settings:

- `dataset_id`
- `ds_prefix`
- `species`
- `ds_config`
- sample IDs
- optional phase-specific settings

Real raw-data locations can either be local relative paths through `data/` or
private paths referenced from `configs_private/`.

## Pipeline phases

| Phase | Main purpose | Representative output |
|---|---|---|
| 01 | alignment / count matrix intake | count matrices or Cell Ranger outputs |
| 02 | QC | cleaned Seurat object |
| 03 | normalization | normalized object |
| 04 | integration | integrated object |
| 05 | clustering | cluster labels and QC plots |
| 06 | annotation | annotated Seurat object |
| 07 | CNV review | CNV labels/scores and annotation feedback |
| 08 | DE/enrichment | pseudobulk DE and enrichment reports |
| 09 | CellChat | communication comparison and summary reports |

The canonical downstream input is expected at:

```text
results/scrna/{ds_prefix}/06_annotate/objects/seurat_annotated_final.rds
```

This object is generated locally and is not tracked by Git.

## Manual gates

Some phases require explicit biological review before downstream steps proceed.
Snakemake gate sentinels live inside `results/` and are not tracked. For example:

```bash
touch results/scrna/<ds_prefix>/06_annotate/.gate_approved
```

Use gates to separate computational completion from biological approval.

## Public release checklist

Before pushing a branch to a public remote, check:

```bash
git status --short
git ls-files | grep -E '(^|/)(data|results|logs|configs_private|Notes|\.snakemake)(/|$)' || true
git ls-files | grep -E '\.(rds|h5ad|h5|bam|fastq.gz|fq.gz|pkl)$' || true
```

Do not push if these commands reveal private data, runtime caches, large objects,
or unpublished handoff notes.
