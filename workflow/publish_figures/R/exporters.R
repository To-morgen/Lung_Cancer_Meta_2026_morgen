# Generic figure-data exporters for Seurat objects.

export_value_or <- function(x, y) if (is.null(x)) y else x

assert_export_string <- function(value, field) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(field, " must be a non-empty string.", call. = FALSE)
  }
  invisible(value)
}

resolve_export_path <- function(project_root, relative_path) {
  assert_export_string(relative_path, "path")
  if (grepl("^(/|[A-Za-z]:[/\\\\])", relative_path)) {
    stop("Export paths must be relative to project_root: ", relative_path, call. = FALSE)
  }

  root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  target <- normalizePath(
    file.path(root, relative_path),
    winslash = "/",
    mustWork = FALSE
  )
  if (!identical(target, root) && !startsWith(target, paste0(root, "/"))) {
    stop("Export path escapes project_root: ", relative_path, call. = FALSE)
  }
  target
}

validate_export_spec <- function(spec) {
  if (!identical(as.integer(export_value_or(spec$schema_version, 0L)), 1L)) {
    stop("Export spec schema_version must be 1.", call. = FALSE)
  }
  assert_export_string(spec$dataset_id, "dataset_id")
  assert_export_string(spec$output_root, "output_root")

  objects <- export_value_or(spec$objects, list())
  if (!is.list(objects) || length(objects) == 0L || is.null(names(objects))) {
    stop("Export spec requires a named objects mapping.", call. = FALSE)
  }
  for (object_id in names(objects)) {
    object_spec <- objects[[object_id]]
    assert_export_string(object_id, "object id")
    assert_export_string(object_spec$path, paste0("objects.", object_id, ".path"))
    assert_export_string(
      object_spec$sample_column,
      paste0("objects.", object_id, ".sample_column")
    )
    assert_export_string(
      object_spec$group_column,
      paste0("objects.", object_id, ".group_column")
    )
  }

  views <- export_value_or(spec$views, list())
  if (!is.list(views) || length(views) == 0L) {
    stop("Export spec requires a non-empty views list.", call. = FALSE)
  }
  view_ids <- vapply(views, function(view) {
    assert_export_string(view$id, "view.id")
    view$id
  }, character(1))
  if (anyDuplicated(view_ids)) {
    stop("View ids must be unique: ", paste(unique(view_ids[duplicated(view_ids)]), collapse = ", "), call. = FALSE)
  }

  allowed_outputs <- c("umap", "dotplot", "composition")
  for (view in views) {
    assert_export_string(view$object, paste0("views.", view$id, ".object"))
    if (!view$object %in% names(objects)) {
      stop("View '", view$id, "' references unknown object '", view$object, "'.", call. = FALSE)
    }
    assert_export_string(
      view$annotation_column,
      paste0("views.", view$id, ".annotation_column")
    )
    assert_export_string(
      view$annotation_level,
      paste0("views.", view$id, ".annotation_level")
    )
    outputs <- unlist(export_value_or(view$outputs, character()), use.names = FALSE)
    if (length(outputs) == 0L || any(!outputs %in% allowed_outputs)) {
      stop(
        "View '", view$id, "' outputs must use: ",
        paste(allowed_outputs, collapse = ", "),
        call. = FALSE
      )
    }
    denominators <- unlist(
      export_value_or(view$composition_denominators, "view"),
      use.names = FALSE
    )
    if (any(!denominators %in% c("view", "object"))) {
      stop("View '", view$id, "' has an invalid composition denominator.", call. = FALSE)
    }
  }

  invisible(spec)
}

load_export_spec <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read an export spec.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Export spec not found: ", path, call. = FALSE)
  }
  spec <- yaml::read_yaml(path)
  validate_export_spec(spec)
  spec
}

validate_object_columns <- function(object, object_spec, view) {
  metadata <- object[[]]
  required <- c(
    object_spec$sample_column,
    object_spec$group_column,
    view$annotation_column
  )
  filters <- export_value_or(view$filters, list())
  if (!is.null(view$subset)) {
    filters <- c(filters, list(view$subset))
  }
  filter_columns <- vapply(filters, function(filter) {
    assert_export_string(filter$column, paste0("views.", view$id, ".filters.column"))
    filter$column
  }, character(1))
  missing <- setdiff(c(required, filter_columns), colnames(metadata))
  if (length(missing) > 0L) {
    stop(
      "View '", view$id, "' is missing object metadata columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

select_view_cells <- function(object, view) {
  metadata <- object[[]]
  keep <- rep(TRUE, nrow(metadata))
  filters <- export_value_or(view$filters, list())
  if (!is.null(view$subset)) {
    filters <- c(filters, list(view$subset))
  }

  for (filter in filters) {
    values <- as.character(metadata[[filter$column]])
    include <- unlist(export_value_or(filter$include, character()), use.names = FALSE)
    exclude <- unlist(export_value_or(filter$exclude, character()), use.names = FALSE)
    if (length(include) > 0L) {
      keep <- keep & values %in% as.character(include)
    }
    if (length(exclude) > 0L) {
      keep <- keep & !values %in% as.character(exclude)
    }
  }

  cells <- rownames(metadata)[keep & !is.na(metadata[[view$annotation_column]])]
  if (length(cells) == 0L) {
    stop("View '", view$id, "' selected zero cells.", call. = FALSE)
  }
  cells
}

resolve_annotation_order <- function(object, cells, view) {
  observed <- unique(as.character(object[[]][cells, view$annotation_column]))
  configured <- unlist(export_value_or(view$annotation_order, character()), use.names = FALSE)
  if (length(configured) == 0L) {
    return(sort(observed))
  }
  missing <- setdiff(observed, configured)
  if (length(missing) > 0L) {
    stop(
      "View '", view$id, "' annotation_order omits: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  as.character(configured[configured %in% observed])
}

flatten_feature_groups <- function(groups) {
  if (is.character(groups)) {
    groups <- list(default = groups)
  }
  if (!is.list(groups) || length(groups) == 0L) {
    stop("Feature definition must be a non-empty list or character vector.", call. = FALSE)
  }

  rows <- list()
  position <- 0L
  for (group_name in names(groups)) {
    features <- unlist(groups[[group_name]], use.names = FALSE)
    for (feature in features) {
      position <- position + 1L
      rows[[position]] <- data.frame(
        feature = as.character(feature),
        feature_group = export_value_or(group_name, "default"),
        feature_order = position,
        stringsAsFactors = FALSE
      )
    }
  }
  feature_table <- do.call(rbind, rows)
  feature_table <- feature_table[!duplicated(feature_table$feature), , drop = FALSE]
  rownames(feature_table) <- NULL
  feature_table
}

load_view_features <- function(view, project_root) {
  if (!is.null(view$feature_source)) {
    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("Package 'yaml' is required to read feature_source.", call. = FALSE)
    }
    source_path <- resolve_export_path(project_root, view$feature_source$path)
    feature_config <- yaml::read_yaml(source_path)
    set_name <- view$feature_source$set
    assert_export_string(set_name, paste0("views.", view$id, ".feature_source.set"))
    groups <- feature_config[[set_name]]
    if (is.null(groups)) {
      stop("Feature set '", set_name, "' not found in ", view$feature_source$path, call. = FALSE)
    }
    return(flatten_feature_groups(groups))
  }
  if (!is.null(view$feature_groups)) {
    return(flatten_feature_groups(view$feature_groups))
  }
  if (!is.null(view$features)) {
    return(flatten_feature_groups(as.character(unlist(view$features, use.names = FALSE))))
  }
  stop("View '", view$id, "' requests dotplot output without features.", call. = FALSE)
}

get_export_matrix <- function(object, assay, layer) {
  assay <- export_value_or(assay, SeuratObject::DefaultAssay(object))
  layer <- export_value_or(layer, "data")
  matrix <- tryCatch(
    SeuratObject::LayerData(object, assay = assay, layer = layer),
    error = function(e) NULL
  )
  if (is.null(matrix)) {
    matrix <- tryCatch(
      SeuratObject::GetAssayData(object, assay = assay, layer = layer),
      error = function(e) SeuratObject::GetAssayData(object, assay = assay, slot = layer)
    )
  }
  matrix
}

atomic_write_csv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".tmp_", tmpdir = dirname(path), fileext = ".csv")
  utils::write.csv(data, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Failed to move temporary output into place: ", path, call. = FALSE)
  }
  invisible(path)
}

export_umap_cells <- function(object, cells, object_spec, view, output_path) {
  reduction <- export_value_or(view$reduction, export_value_or(object_spec$reduction, "umap"))
  embeddings <- SeuratObject::Embeddings(object, reduction = reduction)
  if (ncol(embeddings) < 2L) {
    stop("Reduction '", reduction, "' has fewer than two dimensions.", call. = FALSE)
  }
  missing_cells <- setdiff(cells, rownames(embeddings))
  if (length(missing_cells) > 0L) {
    stop("Reduction '", reduction, "' is missing selected cells.", call. = FALSE)
  }
  metadata <- object[[]][cells, , drop = FALSE]
  output <- data.frame(
    cell_id = cells,
    sample_id = as.character(metadata[[object_spec$sample_column]]),
    group_key = as.character(metadata[[object_spec$group_column]]),
    view_id = view$id,
    annotation_level = view$annotation_level,
    annotation = as.character(metadata[[view$annotation_column]]),
    dim1 = as.numeric(embeddings[cells, 1L]),
    dim2 = as.numeric(embeddings[cells, 2L]),
    reduction = reduction,
    stringsAsFactors = FALSE
  )
  atomic_write_csv(output, output_path)
  output
}

scale_feature_values <- function(values, limits) {
  standard_deviation <- stats::sd(values, na.rm = TRUE)
  scaled <- if (!is.finite(standard_deviation) || standard_deviation == 0) {
    rep(0, length(values))
  } else {
    (values - mean(values, na.rm = TRUE)) / standard_deviation
  }
  pmax(as.numeric(limits[1]), pmin(as.numeric(limits[2]), scaled))
}

export_dotplot_expression <- function(
    object, cells, object_spec, view, project_root, annotation_order, output_path) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required for dotplot export.", call. = FALSE)
  }
  feature_table <- load_view_features(view, project_root)
  assay <- export_value_or(view$assay, export_value_or(object_spec$assay, "RNA"))
  layer <- export_value_or(view$layer, export_value_or(object_spec$layer, "data"))
  expression <- get_export_matrix(object, assay, layer)

  missing_features <- setdiff(feature_table$feature, rownames(expression))
  missing_policy <- export_value_or(view$missing_feature_policy, "error")
  if (length(missing_features) > 0L && identical(missing_policy, "error")) {
    stop(
      "View '", view$id, "' is missing features: ",
      paste(missing_features, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(missing_features) > 0L) {
    warning(
      "View '", view$id, "' skipped missing features: ",
      paste(missing_features, collapse = ", "),
      call. = FALSE
    )
    feature_table <- feature_table[!feature_table$feature %in% missing_features, , drop = FALSE]
  }
  if (nrow(feature_table) == 0L) {
    stop("View '", view$id, "' has no available dotplot features.", call. = FALSE)
  }

  expression <- expression[feature_table$feature, cells, drop = FALSE]
  annotations <- as.character(object[[]][cells, view$annotation_column])
  use_expm1 <- identical(layer, "data") && isTRUE(export_value_or(view$expm1_data, TRUE))
  limits <- unlist(export_value_or(view$scale_limits, c(-2.5, 2.5)), use.names = FALSE)
  if (length(limits) != 2L || limits[1] >= limits[2]) {
    stop("View '", view$id, "' scale_limits must contain two increasing values.", call. = FALSE)
  }

  rows <- list()
  row_index <- 0L
  for (annotation_index in seq_along(annotation_order)) {
    annotation <- annotation_order[[annotation_index]]
    group_cells <- which(annotations == annotation)
    group_matrix <- expression[, group_cells, drop = FALSE]
    average_matrix <- if (use_expm1) expm1(group_matrix) else group_matrix
    averages <- Matrix::rowMeans(average_matrix)
    detected <- Matrix::rowMeans(group_matrix > 0) * 100
    for (feature_index in seq_len(nrow(feature_table))) {
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        view_id = view$id,
        annotation_level = view$annotation_level,
        annotation = annotation,
        annotation_order = annotation_index,
        feature = feature_table$feature[[feature_index]],
        feature_group = feature_table$feature_group[[feature_index]],
        feature_order = feature_table$feature_order[[feature_index]],
        avg_expr = as.numeric(averages[[feature_index]]),
        pct_expr = as.numeric(detected[[feature_index]]),
        n_cells = length(group_cells),
        assay = assay,
        layer = layer,
        avg_expr_method = if (use_expm1) "mean_expm1_data" else "mean_layer_value",
        stringsAsFactors = FALSE
      )
    }
  }
  output <- do.call(rbind, rows)
  output$avg_expr_scaled <- ave(
    output$avg_expr,
    output$feature,
    FUN = function(values) scale_feature_values(values, limits)
  )
  output <- output[, c(
    "view_id", "annotation_level", "annotation", "annotation_order",
    "feature", "feature_group", "feature_order", "avg_expr",
    "avg_expr_scaled", "pct_expr", "n_cells", "assay", "layer",
    "avg_expr_method"
  )]
  atomic_write_csv(output, output_path)
  output
}

sample_group_map <- function(metadata, sample_column, group_column) {
  samples <- as.character(metadata[[sample_column]])
  groups <- as.character(metadata[[group_column]])
  split_groups <- split(groups, samples)
  invalid <- names(split_groups)[vapply(split_groups, function(x) length(unique(x)) != 1L, logical(1))]
  if (length(invalid) > 0L) {
    stop(
      "Samples map to multiple group keys: ",
      paste(invalid, collapse = ", "),
      call. = FALSE
    )
  }
  vapply(split_groups, function(x) unique(x)[1], character(1))
}

export_sample_composition <- function(
    object, cells, object_spec, view, annotation_order, output_path) {
  metadata_all <- object[[]]
  metadata_view <- metadata_all[cells, , drop = FALSE]
  sample_column <- object_spec$sample_column
  group_column <- object_spec$group_column
  sample_groups <- sample_group_map(metadata_all, sample_column, group_column)
  observed_samples <- unique(as.character(metadata_view[[sample_column]]))
  configured_samples <- unlist(export_value_or(view$sample_order, character()), use.names = FALSE)
  sample_order <- if (length(configured_samples) == 0L) {
    sort(observed_samples)
  } else {
    missing <- setdiff(observed_samples, configured_samples)
    if (length(missing) > 0L) {
      stop("View '", view$id, "' sample_order omits: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    as.character(configured_samples[configured_samples %in% observed_samples])
  }

  count_table <- table(
    factor(as.character(metadata_view[[sample_column]]), levels = sample_order),
    factor(as.character(metadata_view[[view$annotation_column]]), levels = annotation_order)
  )
  view_denominators <- table(factor(
    as.character(metadata_view[[sample_column]]),
    levels = sample_order
  ))
  object_denominators <- table(factor(
    as.character(metadata_all[[sample_column]]),
    levels = sample_order
  ))
  denominator_scopes <- unique(as.character(unlist(
    export_value_or(view$composition_denominators, "view"),
    use.names = FALSE
  )))

  rows <- list()
  row_index <- 0L
  for (scope in denominator_scopes) {
    denominators <- if (scope == "view") view_denominators else object_denominators
    for (sample_index in seq_along(sample_order)) {
      sample_id <- sample_order[[sample_index]]
      denominator <- as.integer(denominators[[sample_id]])
      for (annotation_index in seq_along(annotation_order)) {
        annotation <- annotation_order[[annotation_index]]
        n_cells <- as.integer(count_table[sample_id, annotation])
        row_index <- row_index + 1L
        rows[[row_index]] <- data.frame(
          view_id = view$id,
          sample_id = sample_id,
          group_key = unname(sample_groups[[sample_id]]),
          sample_order = sample_index,
          annotation_level = view$annotation_level,
          annotation = annotation,
          annotation_order = annotation_index,
          n_cells = n_cells,
          denominator_cells = denominator,
          denominator_scope = scope,
          proportion = if (denominator > 0L) n_cells / denominator else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  output <- do.call(rbind, rows)
  atomic_write_csv(output, output_path)
  output
}

relative_to_project <- function(path, project_root) {
  root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  sub(paste0("^", root, "/?"), "", normalized)
}

export_figure_data <- function(spec, project_root = ".", output_root = NULL) {
  validate_export_spec(spec)
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Package 'SeuratObject' is required for figure-data export.", call. = FALSE)
  }
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  output_relative <- export_value_or(output_root, spec$output_root)
  output_path <- resolve_export_path(project_root, output_relative)
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

  manifest_rows <- list()
  manifest_index <- 0L
  for (object_id in names(spec$objects)) {
    object_spec <- spec$objects[[object_id]]
    object_views <- Filter(function(view) identical(view$object, object_id), spec$views)
    if (length(object_views) == 0L) next

    source_path <- resolve_export_path(project_root, object_spec$path)
    if (!file.exists(source_path)) {
      stop("Source object not found: ", object_spec$path, call. = FALSE)
    }
    message("Loading object '", object_id, "': ", object_spec$path)
    object <- readRDS(source_path)
    source_info <- file.info(source_path)

    for (view in object_views) {
      validate_object_columns(object, object_spec, view)
      cells <- select_view_cells(object, view)
      annotation_order <- resolve_annotation_order(object, cells, view)
      view_path <- file.path(output_path, view$id)
      dir.create(view_path, recursive = TRUE, showWarnings = FALSE)
      outputs <- as.character(unlist(view$outputs, use.names = FALSE))

      exporters <- list(
        umap = function(path) export_umap_cells(
          object, cells, object_spec, view, path
        ),
        dotplot = function(path) export_dotplot_expression(
          object, cells, object_spec, view, project_root,
          annotation_order, path
        ),
        composition = function(path) export_sample_composition(
          object, cells, object_spec, view, annotation_order, path
        )
      )
      filenames <- c(
        umap = "umap_cells.csv",
        dotplot = "dotplot_expression.csv",
        composition = "sample_composition.csv"
      )

      for (output_type in outputs) {
        target <- file.path(view_path, filenames[[output_type]])
        exported <- exporters[[output_type]](target)
        manifest_index <- manifest_index + 1L
        manifest_rows[[manifest_index]] <- data.frame(
          dataset_id = spec$dataset_id,
          view_id = view$id,
          object_id = object_id,
          source_object = object_spec$path,
          source_size_bytes = as.numeric(source_info$size),
          source_mtime = format(source_info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
          output_type = output_type,
          output_path = relative_to_project(target, project_root),
          n_rows = nrow(exported),
          n_selected_cells = length(cells),
          stringsAsFactors = FALSE
        )
      }
      message("Exported view '", view$id, "' (", length(cells), " cells).")
    }
    rm(object)
    invisible(gc())
  }

  manifest <- do.call(rbind, manifest_rows)
  manifest_path <- file.path(output_path, "export_manifest.csv")
  atomic_write_csv(manifest, manifest_path)
  session_path <- file.path(output_path, "session_info.txt")
  writeLines(capture.output(utils::sessionInfo()), session_path, useBytes = TRUE)
  invisible(manifest)
}
