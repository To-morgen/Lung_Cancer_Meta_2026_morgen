# Dataset Onboarding

This document describes the public-safe path for adding a new lung cancer dataset
to the pipeline.

## 1. Choose dataset class

Use one of the public template classes:

| Class | Use case |
|---|---|
| `public_original` | GEO/SRA/ArrayExpress/TCGA or other public source |
| `internally_generated` | locally generated or collaborator-provided data |
| `vendor_curated` | third-party processed packages |
| `atlas_ready` | already QC-approved standardized objects |

Only public metadata and placeholder paths should be committed.

## 2. Create a dataset metadata config

For public data:

```bash
cp configs/datasets/_template.yaml configs/datasets/<dataset_id>.yaml
```

For internal/private data:

```bash
mkdir -p configs_private/datasets
cp configs/datasets/_TEMPLATE_internally_generated_.yaml \
  configs_private/datasets/<dataset_id>.yaml
```

`configs_private/` is ignored by Git. Keep real vendor delivery paths, sample
manifests, and unpublished sample details there.

## 3. Create a Snakemake entry point

Start from the public example:

```bash
cp -r workflow/datasets/example_lung_2grp workflow/datasets/<dataset_name>
```

Edit `workflow/datasets/<dataset_name>/config.yaml`:

```yaml
dataset_id: "<dataset_id>"
ds_prefix: "<short_output_prefix>"
species: "human"
ds_config: "configs/datasets/<dataset_id>.yaml"
samples:
  - sample_1
  - sample_2
```

Use relative `data/...` paths or references to private configs. Do not commit
absolute local filesystem paths.

## 4. Add optional downstream configs

Annotation, CNV, and DE can each be overridden per dataset:

```yaml
annotation_config: "configs/annotation/celltype_mapping_template.csv"
cnv_config: "configs/annotation/cnv_targets_example.yaml"
de_contrasts_config: "configs/annotation/de_contrasts_example.yaml"
```

For real projects, copy these templates and keep any unpublished biological
interpretation in private notes until it is safe to publish.

## 5. Dry-run before execution

```bash
snakemake -s workflow/datasets/<dataset_name>/Snakefile -n --cores 1
```

A dry-run may report missing input files if example data are not present. That is
expected for templates; it should still parse the Snakefile and list target rules.

## 6. Public release checklist

Before pushing a branch to a public remote:

```bash
git status --short
git ls-files | grep -E '(^|/)(data|results|logs|configs_private|Notes|\.snakemake)(/|$)' || true
git ls-files | grep -E '\.(rds|h5ad|h5|bam|fastq.gz|fq.gz|pkl)$' || true
```

The final tree should contain only reusable code, sanitized configs, and public
documentation.
