# Figure-spec rendering and provenance output.

load_figure_spec <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read a figure spec.", call. = FALSE)
  }
  if (!file.exists(path)) stop("Figure spec not found: ", path, call. = FALSE)
  spec <- yaml::read_yaml(path)
  if (!identical(as.integer(export_value_or(spec$schema_version, 0L)), 1L)) {
    stop("Figure spec schema_version must be 1.", call. = FALSE)
  }
  for (field in c("manuscript_id", "dataset_id")) {
    assert_export_string(spec[[field]], field)
  }
  assert_export_string(spec$output_root, "output_root")
  if (!is.list(spec$display_groups) || length(spec$display_groups) == 0L) {
    stop("Figure spec requires display_groups.", call. = FALSE)
  }
  if (!is.list(spec$figures) || length(spec$figures) == 0L) {
    stop("Figure spec requires figures.", call. = FALSE)
  }
  spec
}

read_panel_source <- function(contract, project_root) {
  if (!identical(contract$source$type, "table")) {
    stop("Panel builders currently require table source contracts.", call. = FALSE)
  }
  path <- resolve_contract_path(project_root, contract$source$path)
  format <- tolower(export_value_or(contract$source$format, tools::file_ext(path)))
  switch(
    format,
    csv = utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
    tsv = utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE),
    txt = utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE),
    stop("Unsupported panel source format: ", format, call. = FALSE)
  )
}

builder_registry <- function() {
  list(
    atlas_umap = build_atlas_umap,
    annotation_dotplot = build_annotation_dotplot,
    sample_composition = build_sample_composition
  )
}

save_panel_plot <- function(plot, path, width, height, dpi = 600) {
  extension <- tolower(tools::file_ext(path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (extension == "pdf") {
    ggplot2::ggsave(
      path, plot = plot, width = width, height = height, units = "in",
      device = grDevices::cairo_pdf, bg = "white"
    )
  } else if (extension == "png" && requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(
      path, plot = plot, width = width, height = height, units = "in",
      device = ragg::agg_png, dpi = dpi, bg = "white"
    )
  } else if (extension %in% c("tif", "tiff") && requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(
      path, plot = plot, width = width, height = height, units = "in",
      device = ragg::agg_tiff, dpi = dpi, bg = "white", compression = "lzw"
    )
  } else {
    ggplot2::ggsave(
      path, plot = plot, width = width, height = height, units = "in",
      dpi = dpi, bg = "white"
    )
  }
  invisible(path)
}

render_figure_spec <- function(
    figure_spec, contract_registry, project_root = ".", output_root = NULL) {
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  output_relative <- export_value_or(output_root, figure_spec$output_root)
  assert_export_string(output_relative, "output_root")
  output_path <- resolve_export_path(project_root, output_relative)
  contracts <- stats::setNames(
    contract_registry$contracts,
    vapply(contract_registry$contracts, function(x) x$id, character(1))
  )
  builders <- builder_registry()
  manifest_rows <- list()
  manifest_index <- 0L

  for (figure in figure_spec$figures) {
    assert_export_string(figure$id, "figure.id")
    for (panel in figure$panels) {
      assert_export_string(panel$id, paste0(figure$id, ".panel.id"))
      assert_export_string(panel$builder, paste0(figure$id, ".", panel$id, ".builder"))
      if (!panel$builder %in% names(builders)) {
        stop("Unknown panel builder: ", panel$builder, call. = FALSE)
      }
      contract_ids <- as.character(unlist(panel$source_contracts, use.names = FALSE))
      if (length(contract_ids) != 1L) {
        stop("Current panel builders require exactly one source contract.", call. = FALSE)
      }
      contract <- contracts[[contract_ids[[1]]]]
      if (is.null(contract)) {
        stop("Unknown source contract: ", contract_ids[[1]], call. = FALSE)
      }
      audit <- audit_source_contract(contract, project_root)
      if (!isTRUE(audit$ready)) {
        stop("Panel source contract is not ready: ", contract$id, " (", audit$reason, ")", call. = FALSE)
      }

      source_data <- read_panel_source(contract, project_root)
      built <- builders[[panel$builder]](
        source_data,
        export_value_or(panel$parameters, list()),
        list(
          display_groups = figure_spec$display_groups,
          dataset_id = figure_spec$dataset_id,
          manuscript_id = figure_spec$manuscript_id
        )
      )
      dimensions <- export_value_or(panel$dimensions, list())
      width <- if (!is.null(dimensions$width_in)) {
        as.numeric(dimensions$width_in)
      } else {
        figure_width(export_value_or(dimensions$layout, "one_half"))
      }
      height <- as.numeric(export_value_or(dimensions$height_in, 3.2))
      formats <- as.character(unlist(export_value_or(panel$output_formats, "pdf"), use.names = FALSE))
      panel_stem <- paste0(panel$id, "_", panel$builder)
      panel_dir <- file.path(output_path, figure$id)
      dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

      source_output <- file.path(panel_dir, paste0(panel_stem, "_source_data.csv"))
      atomic_write_csv(built$source_data, source_output)
      for (format in formats) {
        plot_output <- file.path(panel_dir, paste0(panel_stem, ".", format))
        save_panel_plot(
          built$plot, plot_output, width, height,
          dpi = as.numeric(export_value_or(dimensions$dpi, 600))
        )
      }
      metadata <- list(
        manuscript_id = figure_spec$manuscript_id,
        dataset_id = figure_spec$dataset_id,
        figure_id = figure$id,
        panel_id = panel$id,
        builder = panel$builder,
        source_contracts = contract_ids,
        source_path = contract$source$path,
        dimensions = list(width_in = width, height_in = height),
        output_formats = formats,
        parameters = export_value_or(panel$parameters, list())
      )
      yaml::write_yaml(metadata, file.path(panel_dir, paste0(panel_stem, "_metadata.yaml")))

      manifest_index <- manifest_index + 1L
      manifest_rows[[manifest_index]] <- data.frame(
        figure_id = figure$id,
        panel_id = panel$id,
        builder = panel$builder,
        source_contract = contract$id,
        source_rows = nrow(built$source_data),
        width_in = width,
        height_in = height,
        output_formats = paste(formats, collapse = ";"),
        stringsAsFactors = FALSE
      )
      message("Rendered ", figure$id, panel$id, " with ", panel$builder, ".")
    }
  }
  manifest <- do.call(rbind, manifest_rows)
  atomic_write_csv(manifest, file.path(output_path, "render_manifest.csv"))
  invisible(manifest)
}
