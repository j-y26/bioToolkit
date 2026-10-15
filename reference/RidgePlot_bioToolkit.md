# Ridge plot of genes or per-cell values

Draws the distribution of one or more features in each group of cells as
overlapping density ridges, the first group at the top. Features can be
genes, metadata columns or a mix, each in its own panel with its own
axis. A drop-in alternative to
[`Seurat::RidgePlot()`](https://satijalab.org/seurat/reference/RidgePlot.html)
with a fixed size; it takes the same arguments as
[`VlnPlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/VlnPlot_bioToolkit.md)
where they apply. Needs the ggridges package.

## Usage

``` r
RidgePlot_bioToolkit(
  object,
  features,
  group.by = NULL,
  cells = NULL,
  assay = NULL,
  slot = "data",
  log = FALSE,
  percent = FALSE,
  y.max = NULL,
  zeros = c("include", "exclude", "annotate"),
  overlap = 1.8,
  adjust = 1,
  cols = NULL,
  italic.genes = TRUE,
  ncol = NULL,
  cell.height = 0.3,
  panel.width = 2.5,
  legend.position = "none",
  bg = "transparent"
)
```

## Arguments

- object:

  A Seurat object, or a data frame with one row per cell (e.g. a table
  of QC metrics; row names are the cell names)

- features:

  Features to plot: genes, metadata columns, module scores or embedding
  dimensions (from a data frame: its numeric columns). Gene symbols are
  matched exactly (rat and mouse "Cd3e", human "CD3E"); missing ones are
  left out with a warning

- group.by:

  Column to group cells by: one ridge per group, filled with the group's
  colour. If NULL (default), the active identities of a Seurat object,
  or one group for a data frame

- cells:

  Cells to use, as names or indices. If NULL (default), all

- assay:

  Assay to read genes from. If NULL (default), the default assay; the
  assay and layer used are reported in a message when a feature is a
  gene, and written to the source data

- slot:

  Layer to read genes from, default is "data"

- log:

  Whether to draw values on a log10 axis: TRUE / FALSE for all features,
  or the names of the features to draw so (e.g. \`c("nCount_RNA",
  "nFeature_RNA")\`). Values of 0 or below cannot be drawn and are left
  out, with a message. Default is FALSE

- percent:

  Whether values are fractions (0 to 1) shown as percentages on the
  axis: TRUE / FALSE, or the names of the features (e.g.
  "mito_fraction"). Default is FALSE

- y.max:

  Right end of the value axis (Seurat's name; the value axis of a ridge
  plot is the x-axis): a number, a quantile of each feature's values
  such as "q99", or NULL (default) for the full range. One value for all
  features, or a vector named by feature

- zeros:

  How cells with a value of 0 are handled for genes: "include"
  (default), "exclude" (ridges of the cells above 0 only) or "annotate"
  (as "exclude", with the % of cells above 0 printed at the right end of
  each ridge)

- overlap:

  Height of the tallest ridge, in rows: above 1, ridges overlap the row
  above. Default is 1.8

- adjust:

  Bandwidth multiplier of the densities, default is 1

- cols:

  Colours of the groups: NULL (default) for the colours
  [`DimPlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/DimPlot_bioToolkit.md)
  gives them ("scanpy-20" up to 20 groups), a qualitative palette name,
  or colours named by group (or in level order)

- italic.genes:

  Whether gene names are set in italics, default is TRUE

- ncol:

  Number of columns of panels when \`stack = FALSE\`. If NULL (default),
  as Seurat

- cell.height:

  Distance between ridges (one row) in inches, default is 0.3

- panel.width:

  Length of the value axis in inches, default is 2.5

- legend.position:

  Position of the group legend, default is "none" (the groups are named
  on the axis)

- bg:

  Background colour, default is "transparent"

## Value

A ggplot object for one feature, otherwise a patchwork

## Details

Densities are computed with one bandwidth per feature (Silverman's rule
on all plotted values, times \`adjust\`), so that ridges of a feature
are smoothed alike. Each ridge is scaled to its own height, so the
ridges compare the shapes of the distributions, not the number of cells.
The source data are as in
[`VlnPlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/VlnPlot_bioToolkit.md).

## See also

[`VlnPlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/VlnPlot_bioToolkit.md),
[`save_plot`](https://j-y26.github.io/bioToolkit/reference/save_plot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
RidgePlot_bioToolkit(obj[[]], c("nCount_RNA", "percent.mt"), group.by = "sample",
                     log = "nCount_RNA")
RidgePlot_bioToolkit(obj, "Col1a1", group.by = "cell_type", zeros = "annotate")
} # }
```
