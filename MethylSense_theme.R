# ============================================================================
# MethylSense shared theme and colour configuration
# ============================================================================
# Provides a consistent, publication-ready visual identity for all
# MethylSense plots. Based on hrbrthemes::theme_ipsum() with Arial.
#
# Usage: source this file at the top of any MethylSense script:
#   source(file.path(script_dir, "MethylSense_theme.R"))
#
# Version: 5.7.0
# Date: 2026-03-23
# ============================================================================

if (!requireNamespace("hrbrthemes", quietly = TRUE)) {
  message("[INFO] Installing hrbrthemes for publication-ready plots...")
  install.packages("hrbrthemes", quiet = TRUE)
}
suppressPackageStartupMessages(library(hrbrthemes))

# ============================================================================
# METHYLSENSE COLOUR PALETTE
# ============================================================================

METHYLSENSE_COLOURS <- list(
  # Primary group colours
  Control = "#1B5E20",
  Suspected = "#006064",
  Infected = "#4A148C",
  Case = "#4A148C",

  # Extended palette for multi-group designs
  Extended = c("#1B5E20", "#006064", "#4A148C", "#B71C1C", "#E65100"),

  # Accent colours for specific plot elements
  Hyper = "#D32F2F",
  Hypo = "#1976D2",
  Neutral = "#757575",

  # Ribbon/fill colours (lighter versions)
  Control_light = "#A5D6A7",
  Suspected_light = "#80DEEA",
  Infected_light = "#CE93D8"
)

# Backwards-compatible alias
METHYLSENSE_COLORS <- METHYLSENSE_COLOURS

# ============================================================================
# METHYLSENSE GGPLOT2 THEME
# ============================================================================

#' Publication-ready ggplot2 theme for MethylSense
#'
#' Based on hrbrthemes::theme_ipsum() with Arial font.
#'
#' @param base_size Base font size (default 12)
#' @param grid Which grid lines to show: "XY", "X", "Y", or "" (default "XY")
#' @param axis Which axis lines to show: "xy", "x", "y", or "" (default "xy")
#' @param ... Additional arguments passed to theme_ipsum()
#' @return A ggplot2 theme object
theme_methylsense <- function(base_size = 12, grid = "XY", axis = "xy", ...) {
  hrbrthemes::theme_ipsum(
    base_family = "Arial",
    base_size = base_size,
    plot_title_size = base_size + 4,
    subtitle_size = base_size + 1,
    axis_title_size = base_size + 1,
    caption_size = base_size - 2,
    strip_text_size = base_size + 1,
    grid = grid,
    axis = axis,
    ...
  ) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key.size = ggplot2::unit(1.2, "lines"),
      plot.margin = ggplot2::margin(15, 15, 15, 15),
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )
}
