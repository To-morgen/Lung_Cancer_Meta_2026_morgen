# ============================================================================
# utils_plotting.R — Project-level plotting utilities
# Layer: 1 (project-wide shared themes and helpers)
# ============================================================================

suppressPackageStartupMessages(library(ggplot2))

#' Standard project theme
theme_project <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      strip.text = element_text(face = "bold", size = base_size - 1),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

#' Color palette for groups
group_colors <- function() {
  c("A1" = "#E41A1C", "FL" = "#377EB8", "mc" = "#4DAF4A")
}

#' Before/After fill palette
before_after_colors <- function() {
  c("Before" = "#CCCCCC", "After" = "#4DAF4A")
}

#' Safe ggsave wrapper (creates dir if needed)
safe_ggsave <- function(plot, filename, width = 10, height = 8, dpi = 150, ...) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  
  # Save PDF
  ggsave(sub("\\.[^.]+$", ".pdf", filename), plot, width = width, height = height, ...)
  # Save PNG
  ggsave(sub("\\.[^.]+$", ".png", filename), plot, width = width, height = height, dpi = dpi, ...)
  
  log_msg <- sprintf("[plot] Saved: %s (.pdf + .png)", basename(filename))
  if (exists("log_msg", mode = "function", envir = globalenv())) {
    get("log_msg", envir = globalenv())(log_msg)
  } else {
    cat(log_msg, "\n")
  }
}

cat("[init] scripts/utils/utils_plotting.R loaded\n")
