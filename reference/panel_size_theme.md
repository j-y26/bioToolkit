# Theme giving panels a fixed absolute size

Sets the size of every panel of a plot in inches, so that it is
identical across figures whatever the number of panels and the size or
position of the legend. Add it to any ggplot, then save the plot with
[`save_plot`](https://j-y26.github.io/bioToolkit/reference/save_plot.md),
which sizes the file to fit. Requires ggplot2 \>= 4.0.0.

## Usage

``` r
panel_size_theme(
  panel.width,
  panel.height = panel.width,
  ncol = 1,
  nrow = 1,
  panel.spacing = NULL
)
```

## Arguments

- panel.width:

  Width of a panel in inches: one value for all panel columns, or one
  value per column (e.g. from facets of different numbers of categories)

- panel.height:

  Height of a panel in inches: one value for all panel rows, or one per
  row. Default is \`panel.width\` when it is a single value

- ncol, nrow:

  Number of panel columns and rows, default is 1. Ignored when
  \`panel.width\` (\`panel.height\`) gives one value per column (row)

- panel.spacing:

  Gap between panels in inches. If NULL (default), the theme's own
  spacing is kept

## Value

A ggplot2 theme

## See also

[`save_plot`](https://j-y26.github.io/bioToolkit/reference/save_plot.md),
[`theme_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/theme_bioToolkit.md)

## Examples

``` r
library(ggplot2)
p <- ggplot(mtcars, aes(wt, mpg)) + geom_point() +
    theme_bioToolkit() + panel_size_theme(2.5, 2)
if (FALSE) { # \dontrun{
save_plot(p, "mtcars.pdf")
} # }
```
