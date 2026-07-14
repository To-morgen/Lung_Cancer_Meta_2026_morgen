# Sample-level composition builder with explicit denominator selection.

resolve_group_palette <- function(parameters, display_groups, observed_keys) {
  display_mapping <- unlist(display_groups, use.names = TRUE)
  configured <- unlist(builder_value_or(parameters$group_palette, character()), use.names = TRUE)
  if (length(configured) == 0L) {
    colors <- if (length(observed_keys) == 2L) {
      unname(binary_role_palette())
    } else {
      grDevices::hcl.colors(length(observed_keys), "Dark 3")
    }
    configured <- stats::setNames(colors, observed_keys)
  }
  missing <- setdiff(observed_keys, names(configured))
  if (length(missing) > 0L) {
    stop("group_palette is missing keys: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  stats::setNames(
    unname(configured[observed_keys]),
    unname(display_mapping[observed_keys])
  )
}

build_sample_composition <- function(data, parameters = list(), context = list()) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for sample_composition.", call. = FALSE)
  }
  required <- c(
    "sample_id", "group_key", "annotation", "n_cells",
    "denominator_cells", "denominator_scope", "proportion"
  )
  missing <- setdiff(required, colnames(data))
  if (length(missing) > 0L) {
    stop("sample_composition source is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  denominator_scope <- builder_value_or(parameters$denominator_scope, NULL)
  observed_scopes <- unique(as.character(data$denominator_scope))
  if (is.null(denominator_scope) && length(observed_scopes) != 1L) {
    stop(
      "sample_composition requires denominator_scope when multiple scopes are present: ",
      paste(observed_scopes, collapse = ", "),
      call. = FALSE
    )
  }
  denominator_scope <- builder_value_or(denominator_scope, observed_scopes[[1]])
  data <- data[as.character(data$denominator_scope) == denominator_scope, , drop = FALSE]
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
  data$display_group <- map_display_groups(data$group_key, context$display_groups)
  display_order <- unname(unlist(context$display_groups, use.names = TRUE))
  data$display_group <- factor(data$display_group, levels = display_order)
  observed_keys <- names(unlist(context$display_groups, use.names = TRUE))
  observed_keys <- observed_keys[observed_keys %in% unique(as.character(data$group_key))]
  group_palette <- resolve_group_palette(parameters, context$display_groups, observed_keys)

  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(display_group, proportion, fill = display_group)
  ) +
    ggplot2::geom_point(
      shape = 21,
      size = as.numeric(builder_value_or(parameters$point_size, 2.2)),
      stroke = 0.35,
      colour = "black",
      position = ggplot2::position_jitter(width = 0.06, height = 0, seed = 1)
    ) +
    ggplot2::stat_summary(
      fun = mean,
      geom = "crossbar",
      width = 0.48,
      linewidth = 0.35,
      colour = "black",
      fill = NA
    ) +
    ggplot2::scale_fill_manual(values = group_palette, drop = FALSE) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0.02, 0.1))
    ) +
    ggplot2::labs(x = NULL, y = "Cell proportion", fill = NULL) +
    ggplot2::facet_wrap(
      ~display_annotation,
      ncol = as.integer(builder_value_or(parameters$ncol, 3L)),
      scales = builder_value_or(parameters$facet_scales, "fixed")
    ) +
    theme_cns(
      base_size = as.numeric(builder_value_or(parameters$base_size, 7)),
      legend_position = builder_value_or(parameters$legend_position, "none")
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.border = ggplot2::element_rect(colour = "#BDBDBD", fill = NA, linewidth = 0.3),
      axis.line = ggplot2::element_blank()
    )
  list(plot = plot, source_data = data)
}
