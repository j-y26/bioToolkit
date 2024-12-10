# Some utility helper functions

#' @title Getting default ggplot2 colors
#'
#' @description Since Seurat DimPlot uses ggplot2 default colors when no palette is
#' specified, different categories will result in lack of consistency in
#' colors (ggplot2 ensure the first the last color are the same, but the
#' middle colors are different depending on the number of categories). This
#' function returns the default ggplot2 colors for a given number of categories.
#' This function is directly taken and modified from the following StackOverflow post:
#' https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
#'
#' @param n Number of categories
#'
#' @param h Hue range
#'
#' @param showColor Show the colors, default is FALSE
#'
#' @return A vector of colors (hexadecimal)
#'
#' @importFrom scales show_col
#'
#' @export
#'
#' @references https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
#'
ggplotColours <- function(n = 6, h = c(0, 360) + 15, showColor = FALSE) {
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  colors <- hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)

  # Show the colors, if desired
  if (showColor) {
    scales::show_col(colors)
  } else {
    return(colors)
  }
}


#' @title Mix colors based on weights of each color
#'
#' @description This function mixes colors based on the weights of each color.
#'
#' @param colors A vector of colors (hexadecimal)
#'
#' @param weights A vector of weights (numeric)
#'
#' @return A mixed color (hexadecimal)
#' 
#' @importFrom grDevices col2rgb rgb
#'
#' @export
#'
mix_colors <- function(colors, weights = rep(1, length(colors))) {
  # Check if the length of colors and weights are the same
  if (length(colors) != length(weights)) {
    stop("The length of colors and weights must be the same")
  }

  # Check if colors are hexadecimal
  if (!all(grepl("^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$", colors))) {
    stop("Colors must be hexadecimal")
  }

  # Check if weights are numeric
  if (!all(is.numeric(weights))) {
    stop("Weights must be numeric")
  }

  # Normalize the weights
  weights <- weights / sum(weights)

  # Convert colors to RGB
  colors_rgb <- grDevices::col2rgb(colors)

  # Compute the weighted average of the RGB values
  mixed_color <- colSums(t(colors_rgb) * weights) / sum(weights)

  # Convert the mixed color to hexadecimal
  mixed_color <- grDevices::rgb(mixed_color[1], mixed_color[2], mixed_color[3], maxColorValue = 255)

  return(mixed_color)
}