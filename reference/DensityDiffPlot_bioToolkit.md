# Difference in cell density between conditions on an embedding

Shows where on a 2D embedding (UMAP, tSNE, ...) the cells of a condition
are enriched or depleted relative to a reference condition: the density
of each condition is estimated with a Gaussian kernel, normalized so
that cell numbers do not matter, and their difference is drawn on a
diverging scale. One panel per compared condition, with a fixed size, as
[`DimPlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/DimPlot_bioToolkit.md);
save it with
[`save_plot`](https://j-y26.github.io/bioToolkit/reference/save_plot.md);
the values behind it are returned by
[`export_source_data`](https://j-y26.github.io/bioToolkit/reference/export_source_data.md).

## Usage

``` r
DensityDiffPlot_bioToolkit(
  object,
  split.by,
  reference = NULL,
  compare = NULL,
  sample.by = NULL,
  reduction = "umap",
  dims = c(1, 2),
  cells = NULL,
  full.extent = FALSE,
  measure = c("difference", "log2ratio"),
  draw = c("cells", "raster"),
  density.adjust = 1,
  n.grid = 256,
  mask.quantile = 0.05,
  cols = NULL,
  limits = NULL,
  na.value = "grey88",
  contour = FALSE,
  pt.size = 0.3,
  order = TRUE,
  label = FALSE,
  label.by = NULL,
  label.size = 3.5,
  label.color = "black",
  legend.position = c("right", "left", "top", "bottom", "none"),
  legend.title = NULL,
  ncol = NULL,
  arrow.panels = c("bottom.left", "all", "left", "right", "top", "bottom", "none"),
  panel.width = 4,
  panel.height = panel.width,
  panel.spacing = 0.1,
  raster = NULL,
  raster.threshold = 10000,
  raster.dpi = 1000,
  bg = "transparent",
  seed = 42,
  ...
)
```

## Arguments

- object:

  A Seurat object

- split.by:

  Metadata column of the conditions (e.g. genotype)

- reference:

  Level of \`split.by\` the others are compared with. If NULL (default),
  the first level

- compare:

  Levels of \`split.by\` compared with \`reference\`, one panel each. If
  NULL (default), all other levels

- sample.by:

  Metadata column of the samples. If given, one density is computed per
  sample and the densities of a condition's samples are averaged, so
  that every sample weighs the same and one large sample cannot drive
  the map. If NULL (default), the cells of a condition are pooled

- reduction:

  Name of the reduction, default is "umap"

- dims:

  Two dimensions of the reduction to plot, default is c(1, 2)

- cells:

  Cells to use, as names or indices. If NULL (default), all

- full.extent:

  Whether the axes span all cells of the object rather than only the
  plotted ones, default is FALSE

- measure:

  "difference" (default): the compared density minus the reference
  density, divided by the reference's peak density; or "log2ratio": the
  log2 ratio of the two densities, with 1% of the reference's peak added
  to both

- draw:

  "cells" (default): every cell coloured by the value at its position,
  masked cells in \`na.value\`; or "raster": the values on a grid, over
  the cells in \`na.value\`

- density.adjust:

  Multiplier of the kernel bandwidth (Scott's rule on all plotted
  cells), default is 1. Larger values smooth more

- n.grid:

  Number of grid points along each axis, default is 256

- mask.quantile:

  Regions where both conditions are sparser than this quantile of the
  density at the cells' positions are masked (not coloured), default is
  0.05

- cols:

  Colours of the diverging scale: NULL (default) for "blue-magenta", a
  palette name, or colours from low to high, e.g. \`c(low, mid, high)\`

- limits:

  Cutoff of the colour scale: it runs from -\`limits\` to \`limits\`,
  each end trimmed to the most extreme value drawn on its side, so that
  it gets the full colour. If NULL (default), from the lowest to the
  highest value drawn, centred on 0

- na.value:

  Colour of masked cells, default is "grey88"

- contour:

  Whether to draw contour lines of the density of all plotted cells,
  default is FALSE

- pt.size:

  Size of the cells, default is 0.3

- order:

  Whether cells with the largest absolute values are drawn on top,
  default is TRUE (\`draw = "cells"\`)

- label:

  Whether to label groups of cells at their median position, default is
  FALSE

- label.by:

  Metadata column of the labelled groups. If NULL (default), the active
  identities

- label.size, label.color:

  Size and colour of the labels, default 3.5 and "black"

- legend.position:

  Position of the legend: "right" (default), "left", "top", "bottom" or
  "none"

- legend.title:

  Title of the legend. If NULL (default), the measure and the reference

- ncol:

  Number of columns of panels. If NULL (default), all in one row

- arrow.panels:

  Panels that get the dimension arrows: "bottom.left" (default), "all",
  "left", "right", "top", "bottom" or "none", as in
  [`DimPlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/DimPlot_bioToolkit.md)

- panel.width, panel.height:

  Size of each panel in inches, default 4; \`panel.height\` defaults to
  \`panel.width\`

- panel.spacing:

  Gap between panels in inches, default is 0.1

- raster:

  Whether to rasterize the cells: NULL (default) does so above
  \`raster.threshold\` cells per panel; or TRUE / FALSE

- raster.threshold:

  Number of cells above which \`raster = NULL\` rasterizes, default is
  1e4

- raster.dpi:

  Resolution of rasterized cells, default is 1000

- bg:

  Background colour, default is "transparent"

- seed:

  Random seed of the drawing order of ties and of label placement,
  default is 42

- ...:

  Arguments passed to the dimension arrows, as in
  [`DimPlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/DimPlot_bioToolkit.md)

## Value

A ggplot object

## Details

Cells without a value of \`split.by\` (or \`sample.by\`) are left out,
with a message. All plotted cells set the bandwidth and the axis limits;
cells of conditions that are not compared are drawn but do not enter the
densities.

The map shows where the composition of the embedding shifts between
conditions, not whether the shift is significant: embedding distances
are distorted and cells are not replicates. Use a differential abundance
test on samples (e.g. Milo) for that. With few samples per condition,
consider \`sample.by\` so that one sample cannot dominate.

The source data
([`export_source_data`](https://j-y26.github.io/bioToolkit/reference/export_source_data.md))
hold, for \`draw = "cells"\`, one row per cell and panel: \`cell\`,
\`panel\`, \`x\`, \`y\`, \`condition\`, \`value\`, the normalized
densities \`density_reference\` and \`density_compare\` at the cell, and
whether it is \`masked\`; for \`draw = "raster"\`, the same per grid
point.

## See also

[`DimPlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/DimPlot_bioToolkit.md),
[`FeaturePlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/FeaturePlot_bioToolkit.md),
[`plot_composition_scatter`](https://j-y26.github.io/bioToolkit/reference/plot_composition_scatter.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Where are Mut cells enriched relative to WT?
DensityDiffPlot_bioToolkit(obj, split.by = "genotype", reference = "WT")

# Every sample weighs the same; log2 ratio; group labels
p <- DensityDiffPlot_bioToolkit(obj, split.by = "genotype", reference = "WT",
                                sample.by = "sample", measure = "log2ratio",
                                label = TRUE, label.by = "cell_type")
save_plot(p, "density_difference.pdf")
} # }
```
