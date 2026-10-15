# Generate rank file for GSEA

Ranks the genes of a DE table by their pi-value (see
[`rank_de_genes`](https://j-y26.github.io/bioToolkit/reference/rank_de_genes.md))
and returns, and optionally writes, the two-column table read by the
GSEA desktop application (\`.rnk\`). For fgsea, use
\`setNames(rnk\$rank, rnk\$GeneName)\`.

## Usage

``` r
generate_de_rank_file(
  data,
  file = NULL,
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

- file:

  Path of the tab-separated rank file to write. Default is NULL (nothing
  written)

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

  Deprecated names of the arguments (\`de_results\`, \`save_file\`,
  \`de_method\`, \`pvalue_col\`, \`log2fc_col\`), accepted with a
  warning; renamed in 0.10.0

## Value

A data frame with columns \`GeneName\` (the row names of \`data\`) and
\`rank\` (the pi-value), from highest to lowest; written to \`file\`
when given
