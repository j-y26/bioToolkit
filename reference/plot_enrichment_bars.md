# Diverging bar chart of normalized enrichment scores

Shows the most significant gene sets of one GSEA comparison as bars of
their normalized enrichment score (NES), positive to the right and
negative to the left of 0, ordered by NES. Each gene set's name is
written on the opposite side of 0 from its bar, so that names never
cover bars and no axis labels are needed.

## Usage

``` r
plot_enrichment_bars(
  enrichment,
  comparison = NULL,
  n.up = 10,
  n.down = 10,
  fill = c("direction", "padj"),
  padj.cutoff = 0.05,
  show.ns = FALSE,
  cols = NULL,
  term.replace = NULL,
  wrap.width = NULL,
  cell.height = 0.2,
  panel.width = 4,
  legend.position = "right",
  bg = "transparent",
  ...
)
```

## Arguments

- enrichment:

  Results of
  [`standardize_enrichment`](https://j-y26.github.io/bioToolkit/reference/standardize_enrichment.md),
  or anything it reads, with NES (GSEA results)

- comparison:

  The comparison to draw, when the results hold several. If NULL
  (default), the first

- n.up, n.down:

  Number of gene sets with positive / negative NES to show, the most
  significant first, default is 10 each

- fill:

  Bar colours: "direction" (default; one colour for positive and one for
  negative NES) or "padj" (-log10 adjusted p-value)

- padj.cutoff:

  Adjusted p-value below which a gene set is significant, default is
  0.05

- show.ns:

  Whether non-significant gene sets may be shown (faded, and explained
  in a "Significance" legend), when there are fewer than \`n.up\` /
  \`n.down\` significant ones, default is FALSE

- cols:

  With \`fill = "direction"\`, the colours of positive and negative NES:
  a vector named \`up\` and / or \`down\`, or two colours (up, down).
  NULL (default) takes the ends of the "blue-magenta" palette. With
  \`fill = "padj"\`, a sequential palette name or colours from low to
  high, default "seq-purple"

- term.replace:

  Replacements in gene set names, applied after cleaning: a named
  character vector of regular expressions and their replacements.
  Default is NULL

- wrap.width:

  Names are wrapped onto lines of at most this many characters. If NULL
  (default), as wide as half the panel allows

- cell.height:

  Height of a bar's row in inches, default is 0.2

- panel.width:

  Width of the panel in inches, default is 4

- legend.position:

  Position of the legends, default is "right". With \`fill =
  "direction"\` the colours need no legend (the side of 0 tells the
  direction); the significance legend is drawn when non-significant bars
  are shown

- bg:

  Background colour, default is "transparent"

- ...:

  Arguments to
  [`standardize_enrichment`](https://j-y26.github.io/bioToolkit/reference/standardize_enrichment.md)

## Value

A ggplot object; its data are the plotted gene sets

## See also

[`plot_enrichment_dots`](https://j-y26.github.io/bioToolkit/reference/plot_enrichment_dots.md),
[`standardize_enrichment`](https://j-y26.github.io/bioToolkit/reference/standardize_enrichment.md)

## Examples

``` r
set.seed(1)
fg <- data.frame(pathway = paste0("HALLMARK_SET_", 1:30), pval = runif(30)^4,
                 NES = rnorm(30, sd = 1.5), size = 100)
fg$padj <- p.adjust(fg$pval, "BH")
plot_enrichment_bars(fg, n.up = 5, n.down = 5)

```
