# Split DE results into up- and down-regulated genes

Keeps the genes of a DE table with an adjusted p-value below
\`padj.cutoff\` and an absolute log2 fold change above \`logfc.cutoff\`,
split by direction.

## Usage

``` r
filter_de_genes(
  data,
  logfc.cutoff = 0,
  padj.cutoff = 0.05,
  method = c("auto", "edgeR", "DESeq2", "limma", "seurat"),
  padj.col = NULL,
  logfc.col = NULL,
  ...
)
```

## Arguments

- data:

  A data frame of DE results with a log2 fold change and an adjusted
  p-value column (edgeR, DESeq2, limma or Seurat results, or any table
  with the columns given)

- logfc.cutoff:

  Minimum absolute log2 fold change (strictly greater). Default is 0
  (any change)

- padj.cutoff:

  Maximum adjusted p-value (strictly smaller). Default is 0.05

- method:

  The tool that produced the results: "auto" (default) detects it from
  the column names, as
  [`standardize_de`](https://j-y26.github.io/bioToolkit/reference/standardize_de.md)
  does; or one of "edgeR", "DESeq2", "limma" and "seurat"

- padj.col:

  Column of the (adjusted) p-values used for the ranking. If NULL
  (default), taken from \`method\` (e.g. "FDR" for edgeR, "padj" for
  DESeq2)

- logfc.col:

  Column of the log2 fold changes. If NULL (default), taken from
  \`method\`

- ...:

  Deprecated names of the arguments (\`de_results\`, \`log2fc_min\`,
  \`pvalue_max\`, \`de_method\`, \`pvalue_col\`, \`log2fc_col\`),
  accepted with a warning; renamed in 0.10.0

## Value

A list of two data frames, the rows of \`data\` that pass:

- up_reg:

  up-regulated genes

- down_reg:

  down-regulated genes

## Examples

``` r
de <- data.frame(
    logFC = c(2.1, -1.4, 0.2), PValue = c(1e-8, 1e-4, 0.6),
    FDR = c(3e-8, 2e-4, 0.6), row.names = c("ISG15", "CD3E", "ACTB")
)
filter_de_genes(de, logfc.cutoff = 1)
#> Detected edgeR results (logFC, FDR)
#> $up_reg
#>       logFC PValue   FDR
#> ISG15   2.1  1e-08 3e-08
#> 
#> $down_reg
#>      logFC PValue   FDR
#> CD3E  -1.4  1e-04 2e-04
#> 
```
