# Theme of the reduction plots: the embedding look with fixed panel sizes

The embedding layout of
[`theme_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/theme_bioToolkit.md),
fixed panel sizes from
[`panel_size_theme`](https://j-y26.github.io/bioToolkit/reference/panel_size_theme.md),
and the legend and strip settings of the reduction plots. Text sizes
given as NULL keep the theme's.

## Usage

``` r
embedding_theme(
  panel.width,
  panel.height,
  ncol,
  nrow,
  panel.spacing,
  legend.position,
  legend.text.size,
  strip.text.size,
  strip.text.face,
  bg,
  arrow.scale = 1
)
```

## Arguments

- panel.width, panel.height:

  Size of a single panel in inches

- ncol, nrow:

  Number of panel columns and rows

- panel.spacing:

  Gap between panels in inches

- legend.position:

  Position of the legend

- legend.text.size:

  Font size of the legend text

- strip.text.size, strip.text.face:

  Font size and face of panel labels

- bg:

  Background colour

- arrow.scale:

  The \`arrow.scale\` of the dimension arrows, default is 1

## Value

A ggplot2 theme
