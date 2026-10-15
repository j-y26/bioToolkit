# Grid of pathway enrichment across comparisons

Shows gene sets (rows) across comparisons (columns), e.g. the Hallmark
sets across cell types: colour is the normalized enrichment score (GSEA)
or -log10 adjusted p-value (over-representation), size the -log10
adjusted p-value or the number of genes, and non-significant results are
faded. For the top gene sets of one or a few results, see
[`plot_enrichment_dots`](https://j-y26.github.io/bioToolkit/reference/plot_enrichment_dots.md).

## Usage

``` r
plot_enrichment_grid(
  enrichment,
  terms = NULL,
  n = 5,
  fill = NULL,
  size = NULL,
  padj.cutoff = 0.05,
  significance = c("fade", "grey", "hide"),
  size.max = NULL,
  cluster.terms = TRUE,
  term.replace = NULL,
  wrap.width = 40,
  cols = NULL,
  flip = FALSE,
  cell.width = 0.18,
  cell.height = 0.18,
  legend.position = c("right", "left", "top", "bottom", "none"),
  bg = "transparent",
  ...
)
```

## Arguments

- enrichment:

  Results of
  [`standardize_enrichment`](https://j-y26.github.io/bioToolkit/reference/standardize_enrichment.md),
  or anything it reads (a table, a clusterProfiler result, or a named
  list of them, one per comparison)

- terms:

  Gene sets to show (names after cleaning, or the original
  \`term_id\`s): a character vector, or a named list (groups of rows).
  If NULL (default), the \`n\` most significant of each comparison,
  combined

- n:

  Number of top gene sets per comparison when \`terms\` is NULL, default
  is 5

- fill:

  What the colour shows: "nes" or "padj" (-log10 adjusted p-value). If
  NULL (default), "nes" when the results have NES (GSEA), else "padj"

- size:

  What the dot size shows: "padj" (-log10 adjusted p-value), "count" or
  "gene_ratio". If NULL (default), "padj" when colour shows NES, else
  "count"

- padj.cutoff:

  Adjusted p-value below which a result is significant, default is 0.05

- significance:

  "fade" (default; significant dots opaque with a dark ring, the others
  faded), "grey" (non-significant dots grey) or "hide"

- size.max:

  Value drawn as the largest dot (larger values are drawn at that size).
  If NULL (default), 10 for -log10 adjusted p-values, else the largest
  value

- cluster.terms:

  Whether to order gene sets by hierarchical clustering of their colour
  values across comparisons (within their groups), default is TRUE.
  Otherwise in the order given, or by the first comparison's value

- term.replace:

  Replacements in gene set names, applied after cleaning: a named
  character vector of regular expressions and their replacements, e.g.
  \`c("Extracellular matrix" = "ECM")\`. Default is NULL

- wrap.width:

  Gene set names are wrapped onto lines of at most this many characters,
  default is 40

- cols:

  Colours of the colour scale: NULL (default; "blue-magenta" for NES,
  "seq-purple" for -log10 adjusted p-values), a palette name or colours
  from low to high

- flip:

  Whether to draw comparisons as rows and gene sets as columns, default
  is FALSE

- cell.width, cell.height:

  Width of a column and height of a row in inches, default is 0.18

- legend.position:

  Position of the legends, default is "right"

- bg:

  Background colour, default is "transparent"

- ...:

  Arguments to
  [`standardize_enrichment`](https://j-y26.github.io/bioToolkit/reference/standardize_enrichment.md)
  when \`enrichment\` is not standardized yet (e.g. \`comparison.col\`)

## Value

A ggplot object; its source data
([`export_source_data`](https://j-y26.github.io/bioToolkit/reference/export_source_data.md))
are the plotted rows of the standardized table

## See also

[`plot_enrichment_dots`](https://j-y26.github.io/bioToolkit/reference/plot_enrichment_dots.md),
[`plot_enrichment_bars`](https://j-y26.github.io/bioToolkit/reference/plot_enrichment_bars.md),
[`standardize_enrichment`](https://j-y26.github.io/bioToolkit/reference/standardize_enrichment.md)

## Examples

``` r
set.seed(1)
sets <- paste0("HALLMARK_SET_", 1:12)
gsea <- function() data.frame(pathway = sets, pval = runif(12)^3, NES = rnorm(12, sd = 1.5),
                              size = 100)
res <- lapply(list(A = gsea(), B = gsea(), C = gsea()), function(d) {
    d$padj <- p.adjust(d$pval, "BH")
    d
})
plot_enrichment_grid(res, n = 4)

```
