# Publication Figure Framework

This module provides public-safe infrastructure for building publication
figures from private analysis results. It separates reusable code from
manuscript-specific choices and generated assets.

## Isolation model

| Surface | Location | Git policy |
|---|---|---|
| Generic contracts, builders, themes, schemas, tests | workflow/publish_figures/ | tracked |
| Dataset paths, labels, genes, pathways, panel selection | configs_private/figures/{dataset}/ | ignored |
| Figure-ready tables, panels, composites, manifests | results/publish_figures/{dataset}/{run_id}/ | ignored |
| Rendered manuscript figures | figures/ | ignored until explicit release |

Do not place unpublished findings, real sample identifiers, result values, or
private paths in tracked examples. Synthetic fixtures are the only data allowed
inside this module.

## Contract-first workflow

1. Define private source contracts for the intended panels.
2. Run the audit and classify every source as present, derivable,
   compute_required, decision_required, optional, or blocked.
3. Build deterministic figure-data exporters for derivable contracts.
4. Validate the exported tables before plotting.
5. Render panels from validated tables, not directly from ad hoc objects.
6. Write one source-data table and one provenance record per panel.

## Figure-data export

The generic Seurat exporter reads a private YAML spec and writes one directory
per declared view. A view fixes the object, annotation column, cell filter,
annotation order, and requested output contracts.

    Rscript workflow/publish_figures/export_figure_data.R \
      configs_private/figures/<dataset>/export_spec.yaml \
      "$PROJECT_ROOT"

Real export specs remain ignored. The public schema is documented in
`schemas/export_spec.schema.yaml`; tests generate a synthetic Seurat object at
runtime and do not commit study data.

Supported outputs:

- `umap_cells.csv`: fixed coordinates plus raw sample/group keys and annotation
- `dotplot_expression.csv`: grouped average expression, scaled expression, and
  detection fraction for a private feature panel
- `sample_composition.csv`: zero-complete sample-by-annotation counts and
  proportions with an explicit `denominator_scope` (`view` or `object`)

Filters use declarative `include`/`exclude` values only. Export specs cannot run
arbitrary R expressions. Display labels such as treatment names are added later
by private figure specs, never by the exporter.

## Panel rendering

The initial renderer supports three table-driven builders:

- `atlas_umap`
- `annotation_dotplot`
- `sample_composition`

    Rscript workflow/publish_figures/render_figures.R \
      configs_private/figures/<dataset>/figure_spec.yaml \
      configs_private/figures/<dataset>/source_contracts.yaml \
      "$PROJECT_ROOT"

Each panel writes the rendered files, filtered panel source data, and a metadata
YAML record. Composition panels refuse ambiguous inputs containing multiple
denominator scopes unless the private panel spec selects one explicitly.
Private `annotation_labels` mappings can replace internal annotation keys for
display without changing the raw keys retained in panel source data.

Rendered PDFs or PNGs are not source contracts. For example, a DotPlot is ready
only when its grouped average-expression and detection-fraction table exists.

## Audit command

    Rscript workflow/publish_figures/audit_sources.R \
      configs_private/figures/<dataset>/source_contracts.yaml \
      "$PROJECT_ROOT" \
      results/publish_figures/<dataset>/contract_audit.csv

Add --strict to return a non-zero exit code while required main-figure
contracts are not ready.

## Public example

    Rscript workflow/publish_figures/audit_sources.R \
      workflow/publish_figures/examples/source_contracts.synthetic.yaml \
      . \
      /tmp/publish_figure_contract_audit.csv

    Rscript workflow/publish_figures/tests/test_contracts.R

## Core table interfaces

umap_cells:

    cell_id, sample_id, group_key, annotation_level, annotation, dim1, dim2

dotplot_expression:

    annotation_level, annotation, feature, avg_expr, avg_expr_scaled,
    pct_expr, n_cells, assay, layer

sample_composition:

    sample_id, group_key, annotation_level, annotation, n_cells,
    denominator_cells, denominator_scope, proportion

Internal group keys are mapped to display labels only in the private figure
spec. Generic exporters and builders must not hardcode treatment names.
