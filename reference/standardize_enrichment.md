# Standardize pathway enrichment results

Reads enrichment results from fgsea, clusterProfiler (GSEA or
over-representation) or the GSEA desktop application into one table with
standard column names, for
[`plot_enrichment_dots`](https://j-y26.github.io/bioToolkit/reference/plot_enrichment_dots.md)
and
[`plot_enrichment_bars`](https://j-y26.github.io/bioToolkit/reference/plot_enrichment_bars.md).

## Usage

``` r
standardize_enrichment(
  data,
  method = c("auto", "fgsea", "clusterProfiler", "gsea"),
  comparison.col = NULL,
  clean.terms = TRUE
)
```

## Arguments

- data:

  A data frame of results, a clusterProfiler \`gseaResult\` /
  \`enrichResult\` object, or a named list of either (names =
  comparisons)

- method:

  The tool: "auto" (default) detects it from the columns; or "fgsea",
  "clusterProfiler" (GSEA or over-representation, told apart by their
  columns) or "gsea" (GSEA desktop report, e.g. read with
  `read.delim(..., check.names = FALSE)`; combine its two
  `gsea_report_for_*` files with
  [`rbind()`](https://rdrr.io/r/base/cbind.html) first)

- comparison.col:

  For a single table holding several comparisons, the column naming the
  comparison of each row. Default is NULL

- clean.terms:

  Whether to make MSigDB-style names readable (drop the collection
  prefix such as "HALLMARK\_" or "GOBP\_", underscores to spaces,
  sentence case), default is TRUE. Other names only get a capital first
  letter when their first word is all lower case (GO descriptions). The
  original names are kept in \`term_id\`

## Value

A data frame with columns \`term\`, \`term_id\`, \`comparison\` (when
several), \`nes\`, \`pvalue\`, \`padj\`, \`size\`, \`gene_ratio\`,
\`count\`, \`leading_edge\` and \`type\` ("gsea" or "ora"); columns a
tool does not report are NA

## Details

Columns read:

|  |  |  |  |  |  |
|----|----|----|----|----|----|
| tool | term | nes | padj | size | leading_edge |
| fgsea | pathway | NES | padj | size | leadingEdge |
| clusterProfiler GSEA | Description | NES | p.adjust | setSize | core_enrichment |
| clusterProfiler ORA | Description | none | p.adjust | BgRatio (M of M/N) | geneID |
| GSEA desktop | NAME | NES | FDR q-val | SIZE | none |

For over-representation results, \`gene_ratio\` is k / n from
\`GeneRatio\` ("k/n") and \`count\` is \`Count\`. For the GSEA desktop
application, \`padj\` is its FDR q-value.

## See also

[`plot_enrichment_dots`](https://j-y26.github.io/bioToolkit/reference/plot_enrichment_dots.md),
[`plot_enrichment_bars`](https://j-y26.github.io/bioToolkit/reference/plot_enrichment_bars.md)

## Examples

``` r
fg <- data.frame(
    pathway = c("HALLMARK_HYPOXIA", "HALLMARK_E2F_TARGETS"),
    pval = c(1e-4, 0.2), padj = c(1e-3, 0.3), NES = c(2.1, -0.8), size = c(150, 190)
)
standardize_enrichment(fg)
#>          term              term_id  nes pvalue  padj size gene_ratio count
#> 1     Hypoxia     HALLMARK_HYPOXIA  2.1  1e-04 0.001  150         NA    NA
#> 2 E2F targets HALLMARK_E2F_TARGETS -0.8  2e-01 0.300  190         NA    NA
#>   leading_edge type
#> 1         <NA> gsea
#> 2         <NA> gsea
standardize_enrichment(list("KO vs WT" = fg, "Het vs WT" = fg))
#>          term              term_id comparison  nes pvalue  padj size gene_ratio
#> 1     Hypoxia     HALLMARK_HYPOXIA   KO vs WT  2.1  1e-04 0.001  150         NA
#> 2 E2F targets HALLMARK_E2F_TARGETS   KO vs WT -0.8  2e-01 0.300  190         NA
#> 3     Hypoxia     HALLMARK_HYPOXIA  Het vs WT  2.1  1e-04 0.001  150         NA
#> 4 E2F targets HALLMARK_E2F_TARGETS  Het vs WT -0.8  2e-01 0.300  190         NA
#>   count leading_edge type
#> 1    NA         <NA> gsea
#> 2    NA         <NA> gsea
#> 3    NA         <NA> gsea
#> 4    NA         <NA> gsea
```
