# Publication UMAP builders from validated coordinate tables.

builder_value_or <- function(x, y) if (is.null(x)) y else x

parameter_vector <- function(parameters, field, default = character()) {
  as.character(unlist(builder_value_or(parameters[[field]], default), use.names = FALSE))
}

filter_builder_values <- function(data, column, include = character()) {
  if (length(include) == 0L) return(data)
  data[as.character(data[[column]]) %in% include, , drop = FALSE]
}

annotation_order_from_data <- function(data, configured = character()) {
  observed <- unique(as.character(data$annotation))
  if (length(configured) > 0L) {
    missing <- setdiff(observed, configured)
    if (length(missing) > 0L) {
      stop("annotation_order omits: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    return(configured[configured %in% observed])
  }
  if ("annotation_order" %in% colnames(data)) {
    order_table <- unique(data[, c("annotation", "annotation_order")])
    return(as.character(order_table$annotation[order(order_table$annotation_order)]))
  }
  sort(observed)
}

resolve_named_palette <- function(configured, levels) {
  palette <- unlist(builder_value_or(configured, character()), use.names = TRUE)
  if (length(palette) == 0L) {
    palette <- grDevices::hcl.colors(length(levels), palette = "Dynamic")
    names(palette) <- levels
  }
  missing <- setdiff(levels, names(palette))
  if (length(missing) > 0L) {
    stop("Palette is missing levels: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  palette[levels]
}

map_display_groups <- function(group_keys, display_groups) {
  mapping <- unlist(display_groups, use.names = TRUE)
  missing <- setdiff(unique(as.character(group_keys)), names(mapping))
  if (length(missing) > 0L) {
    stop("display_groups is missing keys: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  unname(mapping[as.character(group_keys)])
}

resolve_annotation_labels <- function(parameters, annotation_order) {
  configured <- unlist(
    builder_value_or(parameters$annotation_labels, character()),
    use.names = TRUE
  )
  labels <- stats::setNames(annotation_order, annotation_order)
  if (length(configured) > 0L) {
    unknown <- setdiff(names(configured), annotation_order)
    if (length(unknown) > 0L) {
      stop("annotation_labels contains unknown keys: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    labels[names(configured)] <- configured
  }
  if (anyDuplicated(unname(labels))) {
    stop("annotation_labels must be unique.", call. = FALSE)
  }
  labels
}

build_atlas_umap <- function(data, parameters = list(), context = list()) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for atlas_umap.", call. = FALSE)
  }
  required <- c("cell_id", "group_key", "annotation", "dim1", "dim2")
  missing <- setdiff(required, colnames(data))
  if (length(missing) > 0L) {
    stop("atlas_umap source is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  data <- filter_builder_values(
    data,
    "annotation",
    parameter_vector(parameters, "annotation_include")
  )
  annotation_order <- annotation_order_from_data(
    data,
    parameter_vector(parameters, "annotation_order")
  )
  data$annotation <- factor(data$annotation, levels = annotation_order)
  annotation_labels <- resolve_annotation_labels(parameters, annotation_order)
  data$display_annotation <- factor(
    unname(annotation_labels[as.character(data$annotation)]),
    levels = unname(annotation_labels)
  )
  palette <- resolve_named_palette(parameters$annotation_palette, annotation_order)
  split_by_group <- isTRUE(builder_value_or(parameters$split_by_group, FALSE))
  if (split_by_group) {
    data$display_group <- map_display_groups(data$group_key, context$display_groups)
    display_order <- unname(unlist(context$display_groups, use.names = TRUE))
    data$display_group <- factor(data$display_group, levels = display_order)
  }

  point_size <- as.numeric(builder_value_or(
    parameters$point_size,
    if (nrow(data) > 30000L) 0.08 else if (nrow(data) > 5000L) 0.14 else 0.28
  ))
  plot <- ggplot2::ggplot(data, ggplot2::aes(dim1, dim2, colour = annotation)) +
    ggplot2::geom_point(
      size = point_size,
      alpha = as.numeric(builder_value_or(parameters$point_alpha, 0.78)),
      stroke = 0
    ) +
    ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2", colour = NULL) +
    theme_cns(
      base_size = as.numeric(builder_value_or(parameters$base_size, 7)),
      legend_position = builder_value_or(parameters$legend_position, "right")
    ) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 0.3),
      legend.key.height = grid::unit(3, "mm")
    )

  if (split_by_group) {
    plot <- plot + ggplot2::facet_wrap(~display_group, nrow = 1, drop = FALSE)
  }

  if (isTRUE(builder_value_or(parameters$label, FALSE))) {
    label_data <- stats::aggregate(cbind(dim1, dim2) ~ annotation, data, stats::median)
    label_data$display_annotation <- unname(
      annotation_labels[as.character(label_data$annotation)]
    )
    counts <- table(data$annotation)
    min_cells <- as.integer(builder_value_or(parameters$min_label_cells, 25L))
    label_data <- label_data[counts[as.character(label_data$annotation)] >= min_cells, , drop = FALSE]
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      plot <- plot + ggrepel::geom_text_repel(
        data = label_data,
        ggplot2::aes(label = display_annotation),
        colour = "black",
        size = as.numeric(builder_value_or(parameters$label_size, 2.1)),
        box.padding = 0.25,
        point.padding = 0.1,
        min.segment.length = 0,
        segment.size = 0.2,
        seed = 1,
        show.legend = FALSE
      )
    } else {
      plot <- plot + ggplot2::geom_text(
        data = label_data,
        ggplot2::aes(label = display_annotation),
        colour = "black",
        size = as.numeric(builder_value_or(parameters$label_size, 2.1)),
        show.legend = FALSE
      )
    }
  }
  list(plot = plot, source_data = data)
}
