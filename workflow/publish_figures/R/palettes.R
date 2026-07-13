# Color-blind-safe semantic palettes. Dataset labels are supplied privately.

binary_role_palette <- function() {
  c(reference = "#4C78A8", intervention = "#D55E00")
}

qualitative_palette <- function(n) {
  colors <- c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
    "#56B4E9", "#F0E442", "#000000", "#8C564B", "#7F7F7F",
    "#17BECF", "#9467BD"
  )
  if (n > length(colors)) {
    stop(
      sprintf("Requested %d colors; the stable palette supports %d.", n, length(colors)),
      call. = FALSE
    )
  }
  colors[seq_len(n)]
}

name_palette <- function(levels, colors = qualitative_palette(length(levels))) {
  levels <- as.character(levels)
  if (anyDuplicated(levels)) {
    stop("Palette levels must be unique.", call. = FALSE)
  }
  if (length(colors) != length(levels)) {
    stop("Palette colors and levels must have equal length.", call. = FALSE)
  }
  stats::setNames(colors, levels)
}
