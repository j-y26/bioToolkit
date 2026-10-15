# Add dimension arrows

Adds two short arrows with the dimension names (e.g. "UMAP-1", "UMAP-2")
to the bottom left corner of an embedding plot and removes its axes, a
compact replacement for full axes. Works on any ggplot, e.g. from
[`Seurat::DimPlot()`](https://satijalab.org/seurat/reference/DimPlot.html);
the package's own embedding plots draw these arrows themselves.

The arrows keep a fixed proportion to the panel they are drawn in: the
sizes below are those of a panel 3.5 inches on its shorter side, and
everything (lengths, line width, label size) is scaled with the panel's
actual size when the plot is drawn. A small panel gets small arrows that
stay clear of the cells, a large one arrows that are not lost in it.
\`arrow.scale\` makes them bigger or smaller at any panel size.

## Usage

``` r
add_dimension_arrows(
  plot,
  reduction = "umap",
  dims = NULL,
  arrow.scale = 1,
  arrow.length = 2,
  x.coord.adjust = 0.25,
  y.coord.adjust = x.coord.adjust,
  arrow.size = 0.25,
  arrow.width = 0.2,
  arrow.line.width = 3,
  arrow.color = "black",
  text.font.size = 12,
  text.dst = 0.3,
  ...
)
```

## Arguments

- plot:

  A ggplot object

- reduction:

  The reduction method used to generate the plot, default is "umap",
  allowed values are "umap", "tsne", "pca", "spatial"

- dims:

  A vector of dimension names, overwrites the default dimension names
  assumed from reduction method

- arrow.scale:

  Size of the arrows relative to the panel, default is 1: arrows 2 cm
  long on a panel of 3.5 inches (about 23% of its shorter side), scaled
  with the panel. Use e.g. 0.8 for smaller or 1.25 for larger arrows

- arrow.length:

  The length of the arrows in cm on a 3.5-inch panel, default is 2

- x.coord.adjust:

  Horizontal adjustment of the arrow origin in cm on a 3.5-inch panel,
  default is 0.25

- y.coord.adjust:

  Vertical adjustment of the arrow origin, default is the same as
  \`x.coord.adjust\`

- arrow.size:

  The length of the arrow head in cm on a 3.5-inch panel, default is
  0.25

- arrow.width:

  The width of the arrow head in cm on a 3.5-inch panel, default is 0.2

- arrow.line.width:

  The width of the arrow line on a 3.5-inch panel (in R's \`lwd\`
  units), default is 3

- arrow.color:

  The color of the arrow and label, default is "black"

- text.font.size:

  The font size of the labels in points on a 3.5-inch panel, default is
  12

- text.dst:

  The distance of the labels from the arrows in cm on a 3.5-inch panel,
  default is 0.3

- ...:

  Deprecated snake_case names of the arguments above (e.g.
  \`arrow_length\`), accepted with a warning; renamed to dot.case in
  0.10.0

## Value

a ggplot object with the dimension arrows added
