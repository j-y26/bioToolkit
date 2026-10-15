# Marker dot plot

Shows, for each gene and group of cells, the percentage of cells
expressing the gene (dot size) and its mean expression (colour), scaled
across groups by default. A drop-in alternative to
[`Seurat::DotPlot()`](https://satijalab.org/seurat/reference/DotPlot.html)
and `SCpubr::do_DotPlot()`, with gene groups, clustering, a coloured
marker for each group, gene names in italics and a fixed size: every row
and column is \`cell.height\` / \`cell.width\` inches, so dots and text
are the same size in every figure. Save it with
[`save_plot`](https://j-y26.github.io/bioToolkit/reference/save_plot.md);
the numbers behind it are returned by
[`export_source_data`](https://j-y26.github.io/bioToolkit/reference/export_source_data.md).

## Usage

``` r
DotPlot_bioToolkit(
  object,
  features,
  group.by = NULL,
  split.by = NULL,
  cells = NULL,
  assay = NULL,
  slot = "data",
  average = c("linear", "log"),
  scale = TRUE,
  expr.threshold = 0,
  pct.min = 0,
  pct.max = 100,
  min.cutoff = NA,
  max.cutoff = NA,
  cluster.features = FALSE,
  cluster.groups = FALSE,
  cols = NULL,
  group.marker = c("dot", "strip", "none"),
  group.cols = NULL,
  feature.group.cols = NULL,
  feature.group.side = c("opposite", "labels"),
  split.by.cols = NULL,
  split.by.side = c("opposite", "labels"),
  group.bar.width = 3.5,
  flip = FALSE,
  italic.genes = TRUE,
  fill.title = NULL,
  size.title = "% expressing",
  cell.width = 0.18,
  cell.height = 0.18,
  legend.position = c("right", "left", "top", "bottom", "none"),
  bg = "transparent"
)
```

## Arguments

- object:

  A Seurat object

- features:

  Genes to plot: a character vector, or a named list of character
  vectors (e.g. markers per cell type, see
  [`top_markers`](https://j-y26.github.io/bioToolkit/reference/top_markers.md)),
  drawn as labelled groups of columns. Genes are matched exactly (rat
  and mouse symbols are e.g. "Cd3e", human "CD3E"); genes not in the
  assay are left out with a warning. A gene in two groups is drawn in
  both

- group.by:

  Metadata column to group cells by: one row per group, marked by the
  group's colour (see \`group.marker\`). If NULL (default), the active
  identities

- split.by:

  Metadata column of a higher-level grouping, drawn like facets: one
  labelled block of rows per level, each with one row per group. E.g.
  \`group.by = "condition", split.by = "cell_type"\` compares conditions
  within each cell type. Default is NULL. (Seurat's \`split.by\` is the
  other way round: it splits each group)

- cells:

  Cells to use, as names or indices. If NULL (default), all

- assay:

  Assay to read expression from. If NULL (default), the default assay of
  the object, as in Seurat and SCpubr. The assay and layer used are
  always reported in a message and written to the source data. Set it
  explicitly when the object holds several assays (e.g. "RNA" and
  "SCT"), which give different values for the same gene

- slot:

  Layer to read, default is "data" (log-normalized values). For "data",
  see \`average\`; other layers are averaged as they are

- average:

  How the mean expression of a group is computed from the log-normalized
  "data" layer:

  "linear"

  :   (default, as
      [`Seurat::DotPlot()`](https://satijalab.org/seurat/reference/DotPlot.html))
      `log1p(mean(expm1(x)))`: the log of the mean normalized
      expression, i.e. a pseudobulk of the group, on the scale of
      Seurat's fold changes

  "log"

  :   (as `SCpubr::do_DotPlot()`) `mean(x)`: the mean of the
      log-normalized values, less influenced by a few cells with very
      high expression

- scale:

  Whether to z-score each gene's mean expression across the plotted
  groups, default is TRUE (as Seurat; SCpubr's \`zscore.data\`). Scaled
  colours compare groups for each gene, not genes with each other. Needs
  at least two groups

- expr.threshold:

  A cell expresses a gene when its value is above this, default is 0

- pct.min:

  Dots of groups with fewer than this percentage of cells expressing the
  gene are not drawn (SCpubr's \`dot.min\`, which defaults to 5),
  default is 0. They still count for scaling and are kept in the source
  data

- pct.max:

  Percentage drawn as the largest dot, default is 100. Lower it for
  genes detected in few cells (e.g. transcription factors); higher
  percentages are drawn at the largest size and the legend says so

- min.cutoff, max.cutoff:

  Ends of the colour scale: a number, a quantile of the plotted values
  such as "q5" / "q95", or NA (default). With \`scale = TRUE\` the scale
  is centred on 0 and each end is the cutoff (NA: +-2.5) or, when closer
  to 0, the most extreme z-score on that side, so that it gets the full
  colour; otherwise 0 and the largest mean. Values beyond the ends get
  the colour of the end

- cluster.features:

  Whether to order genes by hierarchical clustering of their colour
  values (within each group of \`features\`), default is FALSE

- cluster.groups:

  Whether to order the groups by hierarchical clustering of their colour
  values, default is FALSE. With \`split.by\`, the groups are ordered
  the same way in every block, by their values in all blocks

- cols:

  Colours of the expression scale. If NULL (default), "blue-magenta"
  (diverging) with \`scale = TRUE\`, else "viridis_r" (sequential,
  yellow to dark purple). Otherwise a palette name (see
  [`show_palettes`](https://j-y26.github.io/bioToolkit/reference/show_palettes.md))
  or colours from low to high, e.g. \`c("grey90", "darkred")\` or
  \`c(low, mid, high)\`

- group.marker:

  How each group is marked with its colour next to its label: "dot"
  (default; a coloured dot), "strip" (a colour strip) or "none"

- group.cols:

  Colours of the group markers: NULL (default) for the colours
  [`DimPlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/DimPlot_bioToolkit.md)
  gives them ("scanpy-20" up to 20 groups), a qualitative palette name,
  or colours named by group (or in level order)

- feature.group.cols:

  Colours of the bars marking the groups of genes (a named list of
  \`features\`): NULL (default) for black bars, a qualitative palette
  name, or colours named by gene group

- feature.group.side:

  Side of the gene groups: "opposite" (default; opposite the gene names:
  on top, or on the right with \`flip = TRUE\`) or "labels" (beyond the
  gene names)

- split.by.cols:

  Colours of the bars marking the \`split.by\` blocks: NULL (default)
  for black bars, a qualitative palette name, or colours named by level
  of \`split.by\`

- split.by.side:

  Side of the \`split.by\` blocks: "opposite" (default; opposite the
  group names: on the right, or on top with \`flip = TRUE\`) or "labels"
  (beyond the group names)

- group.bar.width:

  Thickness of the bars marking gene groups and \`split.by\` blocks, in
  points, default is 3.5

- flip:

  Whether to draw genes as rows and groups as columns, default is FALSE
  (genes as columns, as Seurat)

- italic.genes:

  Whether gene names are set in italics, default is TRUE. Set FALSE for
  features that are not genes (e.g. proteins of an ADT assay)

- fill.title:

  Title of the colour legend. If NULL (default), "Mean expression", with
  "(z-score)" when scaled

- size.title:

  Title of the size legend, default is "% expressing"

- cell.width, cell.height:

  Width of a column and height of a row in inches, default is 0.18 (13
  pt). Dots scale with them

- legend.position:

  Position of the legends: "right" (default), "left", "top", "bottom" or
  "none"

- bg:

  Background colour, default is "transparent"

## Value

A ggplot object

## Details

Cells without a value of \`group.by\` (or \`split.by\`) are left out,
with a message. Groups (and group-split combinations) with no cells are
not drawn; group colours are still taken from all levels of \`group.by\`
in the whole object, so that they match an embedding of all cells.

The source data
([`export_source_data`](https://j-y26.github.io/bioToolkit/reference/export_source_data.md))
hold one row per group (and split) and gene: the number of cells,
\`pct_expressing\`, \`mean_expression\`, \`mean_expression_scaled\`
(when scaled), whether the dot is \`drawn\`, and the \`assay\`,
\`layer\` and \`average\` used.

## See also

[`top_markers`](https://j-y26.github.io/bioToolkit/reference/top_markers.md),
[`FeaturePlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/FeaturePlot_bioToolkit.md),
[`save_plot`](https://j-y26.github.io/bioToolkit/reference/save_plot.md),
[`export_source_data`](https://j-y26.github.io/bioToolkit/reference/export_source_data.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Genes as a vector, from the SCT assay
DotPlot_bioToolkit(obj, features = c("Cd3e", "Ms4a1", "Cd68"),
                   group.by = "cell_type", assay = "SCT")

# Marker sets per cell type, drawn as labelled groups of columns
markers <- list(
    "T cells" = c("Cd3e", "Cd3d"),
    "B cells" = c("Ms4a1", "Cd79a"),
    "Macrophages" = c("Cd68", "Lyz2")
)
p <- DotPlot_bioToolkit(obj, features = markers, group.by = "cell_type")
save_plot(p, "dotplot.pdf")

# As SCpubr::do_DotPlot(): mean of log values, unscaled, dots from 5\
DotPlot_bioToolkit(obj, features = markers, group.by = "cell_type",
                   average = "log", scale = FALSE, pct.min = 5)

# Top markers from FindAllMarkers()
all_markers <- Seurat::FindAllMarkers(obj, only.pos = TRUE)
DotPlot_bioToolkit(obj, features = top_markers(all_markers, n = 3))
} # }
```
