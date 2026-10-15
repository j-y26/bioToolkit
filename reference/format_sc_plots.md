# Format a set of embedding plots (from Seurat, scCustomize, SCpubr etc.)

Formats a list of ggplot objects or a patchwork object consistently and
arranges them in a grid. By default the axes are removed and replaced by
dimension arrows in the bottom left corner of the bottom left plot;
\`arrow.panels\` puts them on every plot or on one edge of the grid
instead. Arrows a plot already has (e.g. from
[`DimPlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/DimPlot_bioToolkit.md))
are replaced, so that each plot has at most one set.

## Usage

``` r
format_sc_plots(
  plots,
  ncol = NULL,
  nrow = NULL,
  by.row = TRUE,
  arrow.panels = c("bottom.left", "all", "left", "right", "top", "bottom", "none"),
  remove.axes = TRUE,
  remove.legend = FALSE,
  remove.grid = FALSE,
  restore.legend = FALSE,
  legend.position = NULL,
  legend.title = NULL,
  legend.title.size = 12,
  legend.title.bold = TRUE,
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

- plots:

  A list of ggplot objects or a patchwork object

- ncol:

  Number of columns in the plot grid, default is NULL. Prioritized if
  nrow is also set.

- nrow:

  Number of rows in the plot grid, default is NULL

- by.row:

  Whether to fill the plot grid by row, default is TRUE, consistent with
  patchwork's default behavior. If FALSE, the plots will be filled by
  column.

- arrow.panels:

  Plots that get the dimension arrows: "bottom.left" (default; the first
  plot of the bottom row), "all", "left" / "right" (the first / last
  plot of each row), "top" / "bottom" (the top / lowest plot of each
  column) or "none". Edges are those of the plots drawn: in a grid whose
  last row is incomplete, a plot with nothing below it is on the bottom
  edge. The grid is given by \`ncol\` / \`nrow\`, the patchwork's
  layout, or patchwork's default. Unless "none", the axes of all plots
  are removed. (Replaces \`add.arrow\` and \`bottom.left.only\`,
  accepted with a deprecation warning)

- remove.axes:

  Whether to remove the axes, default is TRUE

- remove.legend:

  Whether to remove the legend, default is FALSE

- remove.grid:

  Whether to remove the grid, default is FALSE

- restore.legend:

  Whether to restore the legend, default is FALSE

- legend.position:

  The position of the legend, default is NULL (meaning no change)

- legend.title:

  The title of the legend, default is NULL

- legend.title.size:

  The size of the legend title, default is 12

- legend.title.bold:

  Whether to bold the legend title, default is TRUE

- reduction:

  The reduction method used to generate the plots, default is "umap",
  allowed values are "umap", "tsne", "pca", "spatial". Used only when
  arrows are drawn

- dims:

  A vector of dimension names, overwrites the default dimension names.
  Used only when arrows are drawn

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

  Deprecated arguments, accepted with a warning: \`add.arrow\` and
  \`bottom.left.only\` (replaced by \`arrow.panels\`), and the
  snake_case names renamed to dot.case in 0.10.0 (e.g. \`byrow\`,
  \`legend_position\`, \`arrow_length\`)

## Value

a patchwork object with the formatted plots
