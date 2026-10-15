# Getting default ggplot2 colors

Since Seurat DimPlot uses ggplot2 default colors when no palette is
specified, different categories will result in lack of consistency in
colors (ggplot2 ensure the first the last color are the same, but the
middle colors are different depending on the number of categories). This
function returns the default ggplot2 colors for a given number of
categories. This function is directly taken and modified from the
following StackOverflow post:
https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette

## Usage

``` r
ggplot_colors(n = 6, h = c(0, 360) + 15, show.colors = FALSE)
```

## Arguments

- n:

  Number of categories

- h:

  Hue range

- show.colors:

  Whether to draw the colors (with
  [`scales::show_col()`](https://scales.r-lib.org/reference/show_col.html))
  instead of returning them, default is FALSE

## Value

A vector of colors (hexadecimal)

## References

https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
