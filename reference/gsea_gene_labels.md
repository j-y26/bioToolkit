# Names of a few genes of a GSEA plot, with an arrow to their ticks

The names are one comma-separated line (wrapped when wider than 80% of
the panel, measured on a PDF device with 15% to spare), a little inside
the panel from the end of the list the genes are at, and a small arrow
runs from the names to the edge of the panel next to the gene ticks, at
the middle rank of the named genes.

## Usage

``` r
gsea_gene_labels(
  genes,
  target,
  n,
  positive,
  y_rng,
  to,
  panel.width,
  panel.height,
  italic,
  sz
)
```

## Arguments

- genes:

  Genes named, in order

- target:

  Rank the arrow points to

- n:

  Number of genes in the list

- positive:

  Whether the score is positive (genes at the start of the list)

- y_rng:

  Limits of the panel's y-axis

- to:

  "top" (the arrow points up, from the metric panel) or "bottom" (down,
  from the score panel)

- panel.width, panel.height:

  Size of the panel in inches

- italic:

  Whether the names are set in italics

- sz:

  Sizes from
  [`theme_sizes()`](https://j-y26.github.io/bioToolkit/reference/theme_sizes.md)

## Value

A list of ggplot2 layers
