# Calculate ranks of genes for differential expression analysis

Ranks the genes of a DE table by a signed significance score, the
pi-value: \$\$pi = -log10(padj) \* sign(logFC)\$\$ so that significantly
up-regulated genes come first and significantly down-regulated genes
last; a ranking for GSEA (see
[`generate_de_rank_file`](https://j-y26.github.io/bioToolkit/reference/generate_de_rank_file.md)).

## Usage

``` r
rank_de_genes(
  data,
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

  Deprecated names of the arguments (\`de_results\`, \`de_method\`,
  \`pvalue_col\`, \`log2fc_col\`), accepted with a warning; renamed in
  0.10.0

## Value

The input data frame, ordered by rank, with two added columns:

- piValue:

  The pi-value of the gene

- rank:

  The rank of the gene based on the pi-value

## Details

Ties are common among genes with the same adjusted p-value (e.g. many
non-significant genes at padj = 1); they are broken by row order
(\`ties.method = "first"\`).

## Examples

``` r
de <- data.frame(
    logFC = c(2.1, -1.4, 0.2), PValue = c(1e-8, 1e-4, 0.6),
    FDR = c(3e-8, 2e-4, 0.6), row.names = c("ISG15", "CD3E", "ACTB")
)
rank_de_genes(de)
#> Detected edgeR results (logFC, FDR)
#>       logFC PValue   FDR    piValue rank
#> ISG15   2.1  1e-08 3e-08  7.5228787    1
#> ACTB    0.2  6e-01 6e-01  0.2218487    2
#> CD3E   -1.4  1e-04 2e-04 -3.6989700    3
```
