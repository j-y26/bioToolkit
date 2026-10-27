# Volcano plot of differential expression results

Plots the log2 fold change of each feature against its -log10 adjusted
p-value, colouring features that are significantly up- or
down-regulated, with optional labels and counts of DE features. Reads
edgeR, DESeq2, limma and Seurat results directly (see
[`standardize_de`](https://j-y26.github.io/bioToolkit/reference/standardize_de.md)).

## Usage

``` r
plot_volcano(
  data,
  method = c("auto", "edgeR", "DESeq2", "limma", "seurat"),
  feature.col = NULL,
  logfc.col = NULL,
  p.col = NULL,
  padj.col = NULL,
  padj.cutoff = 0.05,
  logfc.cutoff = 0,
  label = c("top", "selected", "none"),
  label.n = 10,
  label.features = NULL,
  label.rank = c("padj", "logfc", "product"),
  label.size = NULL,
  italic.genes = TRUE,
  annotate.counts = c("corner", "arrow", "none"),
  symmetric = TRUE,
  cols = NULL,
  up.label = "Up-regulated",
  down.label = "Down-regulated",
  ns.label = "Not significant",
  pt.size = 1,
  alpha = 0.6,
  raster = NULL,
  raster.dpi = 1000,
  panel.width = 3,
  panel.height = panel.width,
  legend.position = NULL,
  bg = "transparent",
  seed = 42
)
```

## Arguments

- data:

  DE results of one comparison: a data frame (or a DESeq2
  \`DESeqResults\` / edgeR \`TopTags\`), with feature names as row names
  or in \`feature.col\`

- method:

  The tool that produced the results: "auto" (default) detects it from
  the column names; or one of "edgeR", "DESeq2", "limma" and "seurat"

- feature.col:

  Column of feature (gene) names. If NULL (default), the row names (or
  the \`gene\` column of Seurat results)

- logfc.col:

  Column of the log2 fold change. If NULL (default), taken from
  \`method\`. To plot a shrunken log fold change (DESeq2 `lfcShrink()`,
  e.g. with apeglm), which is less noisy for lowly expressed genes, give
  its column here

- p.col:

  Column of the raw p-value. If NULL (default), taken from \`method\`.
  Only used to compute adjusted p-values (Benjamini-Hochberg) when the
  table has none

- padj.col:

  Column of the adjusted p-value. If NULL (default), taken from
  \`method\`. It is both plotted on the y-axis and used to decide
  significance

- padj.cutoff:

  Adjusted p-value below which a feature is significant, default is
  0.05. Drawn as a dashed horizontal line

- logfc.cutoff:

  Absolute log2 fold change a significant feature must exceed, default
  is 0 (any change). When above 0, drawn as dashed vertical lines

- label:

  Which features to label:

  "top"

  :   (default) the \`label.n\` most significant features of each
      direction, ranked by \`label.rank\`, plus any \`label.features\`

  "selected"

  :   only \`label.features\`, whether or not they are significant

  "none"

  :   no labels

- label.n:

  Number of features to label per direction with \`label = "top"\`,
  default is 10. Every requested label is drawn; when long labels (e.g.
  Ensembl IDs) crowd the same corner they can overlap, so lower
  \`label.n\` or widen \`panel.width\`

- label.features:

  Features to label. Always labelled with \`label = "selected"\`, and in
  addition to the top features with \`label = "top"\`. Matched exactly,
  so the case of gene symbols matters (e.g. "CD3E" in human, "Cd3e" in
  mouse and rat)

- label.rank:

  How the top features are ranked with \`label = "top"\`: "padj"
  (default; smallest adjusted p-value, ties broken by the largest
  absolute log2 fold change), "logfc" (largest absolute log2 fold
  change) or "product" (-log10 adjusted p-value x absolute log2 fold
  change, which favours features that are both)

- label.size:

  Size of the labels in mm. If NULL (default), the small small text size
  of the theme (9 pt)

- italic.genes:

  Whether labels are gene names, set in italics. Default is TRUE; set
  FALSE for features that are not genes (e.g. peaks or pathways). Labels
  are bold either way

- annotate.counts:

  How to show the numbers of up- and down-regulated features: "corner"
  (default; in the bottom corners), "arrow" (arrows naming the direction
  at the top of the panel, with the counts below them) or "none"

- symmetric:

  Whether the x-axis is symmetric around 0, default is TRUE

- cols:

  Colours of the up-regulated, down-regulated and not significant
  features: a vector named \`up\`, \`down\` and \`ns\` (any of them; the
  others keep their default), or three colours in that order. NULL
  (default) is \`c(up = "#C91F28", down = "#174EA6", ns = "grey75")\`

- up.label, down.label, ns.label:

  Names of the three classes in the legend (and next to the arrows with
  \`annotate.counts = "arrow"\`). Default "Up-regulated",
  "Down-regulated" and "Not significant"

- pt.size:

  Point size, default is 1

- alpha:

  Point opacity, default is 0.6

- raster:

  Whether to rasterize the points while keeping labels, lines and text
  as vectors: NULL (default) rasterizes when more than 5,000 features
  are drawn, or TRUE / FALSE

- raster.dpi:

  Resolution of the rasterized points, default is 1000. Match the
  \`dpi\` used when saving raster (e.g. png) output

- panel.width:

  Width of the panel in inches, default is 3

- panel.height:

  Height of the panel in inches. Default is \`panel.width\`, a square
  panel, so that setting \`panel.width\` alone resizes the plot

- legend.position:

  Position of the legend. If NULL (default), "right", or "none" with
  \`annotate.counts = "arrow"\` (the arrows already name the
  directions). A legend at the "top" is usually wider than the panel,
  which widens the saved figure

- bg:

  Background colour, default is "transparent"

- seed:

  Random seed for the placement of labels, default is 42

## Value

A ggplot object

## Details

A feature is up-regulated when its adjusted p-value is below
\`padj.cutoff\` and its log2 fold change is above \`logfc.cutoff\`,
down-regulated when below \`-logfc.cutoff\`, and not significant
otherwise. Features without a log fold change or an adjusted p-value are
not drawn.

The plot has a fixed panel size; save it with
[`save_plot`](https://j-y26.github.io/bioToolkit/reference/save_plot.md),
which sizes the file to fit. Its data (one row per drawn feature, with
columns \`direction\`, \`neg_log10_padj\` and \`labelled\` added to the
standardized DE table) can be written with
[`export_source_data`](https://j-y26.github.io/bioToolkit/reference/export_source_data.md).

## See also

[`standardize_de`](https://j-y26.github.io/bioToolkit/reference/standardize_de.md),
[`save_plot`](https://j-y26.github.io/bioToolkit/reference/save_plot.md),
[`export_source_data`](https://j-y26.github.io/bioToolkit/reference/export_source_data.md)

## Examples

``` r
set.seed(1)
de <- data.frame(
    logFC = rnorm(2000, sd = 1.5),
    logCPM = runif(2000, 0, 10),
    PValue = runif(2000)^4,
    row.names = paste0("Gene", 1:2000)
)
de$FDR <- p.adjust(de$PValue, "BH")
plot_volcano(de)
#> Detected edgeR results (logFC, PValue, FDR)

plot_volcano(de, label = "selected", label.features = c("Gene1", "Gene2"))
#> Detected edgeR results (logFC, PValue, FDR)

plot_volcano(de, annotate.counts = "arrow", logfc.cutoff = 1)
#> Detected edgeR results (logFC, PValue, FDR)

if (FALSE) { # \dontrun{
save_plot(plot_volcano(de), "volcano.pdf")
} # }
```
