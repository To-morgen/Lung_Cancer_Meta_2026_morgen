# Marker DotPlot builder from grouped expression tables.

build_annotation_dotplot <- function(data, parameters = list(), context = list()) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for annotation_dotplot.", call. = FALSE)
  }
  required <- c("annotation", "feature", "avg_expr_scaled", "pct_expr", "n_cells")
  missing <- setdiff(required, colnames(data))
  if (length(missing) > 0L) {
    stop("annotation_dotplot source is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  annotation_include <- parameter_vector(parameters, "annotation_include")
  feature_include <- parameter_vector(parameters, "feature_include")
  data <- filter_builder_values(data, "annotation", annotation_include)
  data <- filter_builder_values(data, "feature", feature_include)
  feature_group_include <- parameter_vector(parameters, "feature_group_include")
  if (length(feature_group_include) > 0L && "feature_group" %in% colnames(data)) {
    data <- filter_builder_values(data, "feature_group", feature_group_include)
  }

  annotation_order <- annotation_order_from_data(
    data,
    parameter_vector(parameters, "annotation_order")
  )
  feature_order <- if (length(feature_include) > 0L) {
    feature_include[feature_include %in% unique(as.character(data$feature))]
  } else if ("feature_order" %in% colnames(data)) {
    feature_table <- unique(data[, c("feature", "feature_order")])
    as.character(feature_table$feature[order(feature_table$feature_order)])
  } else {
    unique(as.character(data$feature))
  }
  annotation_labels <- resolve_annotation_labels(parameters, annotation_order)
  data$annotation <- factor(data$annotation, levels = annotation_order)
  data$display_annotation <- factor(
    unname(annotation_labels[as.character(data$annotation)]),
    levels = rev(unname(annotation_labels))
  )
  data$feature <- factor(data$feature, levels = feature_order)
  if ("feature_group" %in% colnames(data)) {
    feature_table <- unique(data[, c("feature", "feature_group", "feature_order")])
    group_order <- unique(as.character(feature_table$feature_group[order(feature_table$feature_order)]))
    data$feature_group <- factor(data$feature_group, levels = group_order)
  }

  color_limits <- as.numeric(unlist(builder_value_or(parameters$color_limits, c(-2.5, 2.5))))
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(feature, display_annotation, size = pct_expr, colour = avg_expr_scaled)
  ) +
    ggplot2::geom_point(alpha = 0.95) +
    ggplot2::scale_size_continuous(
      name = "% expressing",
      range = as.numeric(unlist(builder_value_or(parameters$size_range, c(0.2, 3.4)))),
      limits = c(0, 100),
      breaks = c(25, 50, 75, 100)
    ) +
    ggplot2::scale_colour_gradient2(
      name = "Scaled\nexpression",
      low = builder_value_or(parameters$low_color, "#2166AC"),
      mid = builder_value_or(parameters$mid_color, "#F7F7F7"),
      high = builder_value_or(parameters$high_color, "#B2182B"),
      midpoint = 0,
      limits = color_limits,
      oob = scales::squish
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_cns(
      base_size = as.numeric(builder_value_or(parameters$base_size, 7)),
      legend_position = builder_value_or(parameters$legend_position, "right")
    ) +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1),
      panel.grid.major = ggplot2::element_line(colour = "#E5E5E5", linewidth = 0.25),
      panel.grid.minor = ggplot2::element_blank(),
      strip.placement = "outside"
    )

  if ("feature_group" %in% colnames(data) &&
      isTRUE(builder_value_or(parameters$facet_feature_groups, TRUE))) {
    plot <- plot + ggplot2::facet_grid(
      cols = ggplot2::vars(feature_group),
      scales = "free_x",
      space = "free_x",
      switch = "x"
    )
  }
  list(plot = plot, source_data = data)
}
