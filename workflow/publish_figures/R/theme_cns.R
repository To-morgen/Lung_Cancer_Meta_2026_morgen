# Compact publication theme and physical-size registry.

figure_width <- function(layout = c("single", "one_half", "double")) {
  layout <- match.arg(layout)
  unname(c(single = 3.35, one_half = 5.2, double = 7.2)[layout])
}

theme_cns <- function(
    base_size = 7,
    base_family = "sans",
    legend_position = "right") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for theme_cns().", call. = FALSE)
  }

  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 0.5, colour = "black"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base_size, face = "bold"),
      legend.position = legend_position,
      legend.title = ggplot2::element_text(size = base_size),
      legend.text = ggplot2::element_text(size = base_size - 0.5),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(3.5, "mm"),
      panel.spacing = grid::unit(1.5, "mm"),
      plot.margin = ggplot2::margin(2, 2, 2, 2, unit = "mm")
    )
}
