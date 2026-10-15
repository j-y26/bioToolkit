# Columns of a DE table used by the DE table helpers

Finds the log2 fold change and adjusted p-value columns of a DE table
for
[`rank_de_genes()`](https://j-y26.github.io/bioToolkit/reference/rank_de_genes.md),
[`generate_de_rank_file()`](https://j-y26.github.io/bioToolkit/reference/generate_de_rank_file.md)
and
[`filter_de_genes()`](https://j-y26.github.io/bioToolkit/reference/filter_de_genes.md),
from the same presets as
[`standardize_de`](https://j-y26.github.io/bioToolkit/reference/standardize_de.md):
columns given explicitly win; otherwise they come from \`method\`,
detected from the column names when "auto".

## Usage

``` r
de_helper_columns(data, method, logfc.col, padj.col)
```

## Arguments

- data:

  A DE table

- method:

  "auto", "edgeR", "DESeq2", "limma" or "seurat" (case is ignored, as
  the old \`de_method\` argument did)

- logfc.col, padj.col:

  Columns given by the user, or NULL

## Value

A named character vector with \`logfc\` and \`padj\`
