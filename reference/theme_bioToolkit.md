# The bioToolkit ggplot2 theme

The theme shared by all plots of the package: 11 pt text (9 pt tick
labels and legend text, 12 pt titles), 1 pt axis lines and 0.5 pt other
lines. Three layouts are provided:

- "axis":

  Classic axes: axis lines and black ticks, no grid. For scatter,
  volcano, MA, bar and PCA plots.

- "grid":

  The axes of "axis", plus light grid lines, no axis titles and x labels
  rotated by 45 degrees. For dot plots and heatmaps, whose rows and
  columns are categories.

- "embedding":

  No axes at all, with extra room at the bottom left for the dimension
  arrows. For UMAP / tSNE plots.

Panel sizes are not part of the theme; add
[`panel_size_theme`](https://j-y26.github.io/bioToolkit/reference/panel_size_theme.md)
to fix them.

## Usage

``` r
theme_bioToolkit(
  type = c("axis", "grid", "embedding"),
  base.size = NULL,
  base.family = NULL,
  bg = "transparent"
)
```

## Arguments

- type:

  Layout of the theme, one of "axis" (default), "grid" or "embedding"

- base.size:

  Base font size in points. If NULL (default), 11; otherwise the other
  text sizes are scaled with it

- base.family:

  Font family. If NULL (default), the device default (Helvetica in PDF,
  which journals accept as Arial's equivalent)

- bg:

  Background colour, default is "transparent" (NA also means
  transparent), so that figures can be placed on any background. Use
  e.g. "white" for an opaque background

## Value

A ggplot2 theme

## See also

[`panel_size_theme`](https://j-y26.github.io/bioToolkit/reference/panel_size_theme.md),
[`save_plot`](https://j-y26.github.io/bioToolkit/reference/save_plot.md),
[`italic_labels`](https://j-y26.github.io/bioToolkit/reference/italic_labels.md)

## Examples

``` r
library(ggplot2)
p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
p + theme_bioToolkit()

p + theme_bioToolkit(base.size = 8)

```
