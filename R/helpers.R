# Some utility helper functions

#' Getting default ggplot2 colors
#' 
#' Since Seurat DimPlot uses ggplot2 default colors when no palette is
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