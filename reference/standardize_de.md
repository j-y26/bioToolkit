# Standardize differential expression results

Reads DE results from edgeR, DESeq2, limma or Seurat (or any table, with
the columns given explicitly) into one table with standard column names,
so that every DE plot of the package reads them the same way.

## Usage

``` r
standardize_de(
  data,
  method = c("auto", "edgeR", "DESeq2", "limma", "seurat"),
  feature.col = NULL,
  logfc.col = NULL,
  p.col = NULL,
  padj.col = NULL,
  mean.col = NULL,
  comparison.col = NULL,
  groups = NULL
)
```

## Arguments

- data:

  A data frame of DE results (or an object coercible to one, such as a
  DESeq2 \`DESeqResults\` or edgeR \`TopTags\`); or a named list of
  them, one per comparison, the names becoming the \`comparison\` column

- method:

  The tool that produced the results: "auto" (default) detects it from
  the column names and says which was picked; or one of "edgeR",
  "DESeq2", "limma" and "seurat"

- feature.col:

  Column of feature (gene) names. If NULL (default), the \`gene\` column
  of Seurat results when there is one (as from `FindAllMarkers()`), else
  the row names

- logfc.col, p.col, padj.col, mean.col:

  Columns of the log2 fold change, raw p-value, adjusted p-value and
  mean expression. Each NULL (default) is taken from \`method\`. Give
  e.g. the column of a shrunken log fold change (DESeq2 `lfcShrink()`,
  apeglm) as \`logfc.col\`

- comparison.col:

  For a single table holding several comparisons, the column naming the
  comparison of each row (e.g. "cluster" for `FindAllMarkers()`).
  Default is NULL

- groups:

  Optional second level of comparisons, used for brackets in
  multi-comparison plots: a data frame with \`comparison\` and \`group\`
  columns, or a character vector of groups named by comparison

## Value

A data frame with the columns \`feature\`, \`logfc\`, \`pvalue\`,
\`padj\`, \`mean_expr\` (when available), \`comparison\` and \`group\`
(when given), then the other columns of the input. Attributes \`method\`
(one per comparison) and \`mean.scale\` ("log2" or "linear", NA when
unknown) record how it was read

## Details

Columns read for each method:

|        |                |         |           |                   |
|--------|----------------|---------|-----------|-------------------|
| method | logfc          | pvalue  | padj      | mean_expr         |
| edgeR  | logFC          | PValue  | FDR       | logCPM (log2)     |
| DESeq2 | log2FoldChange | pvalue  | padj      | baseMean (linear) |
| limma  | logFC          | P.Value | adj.P.Val | AveExpr (log2)    |
| seurat | avg_log2FC     | p_val   | p_val_adj | none              |

When there is no adjusted p-value column, adjusted p-values are computed
with the Benjamini-Hochberg method, per comparison, with a message.

p-values of exactly 0, which come from numerical underflow, are replaced
by the smallest non-zero p-value of the table times 0.1, with a message,
so that they are drawn just beyond the other points rather than at an
arbitrary 1e-300. NA p-values (e.g. from DESeq2's independent filtering)
are kept and reported; plots treat them as not significant.

Seurat versions before 4 reported natural-log fold changes
(\`avg_logFC\`); these are not read automatically, since mixing them
with log2 values would be wrong. Convert them (divide by log(2)) and
pass the column as \`logfc.col\`.

## Examples

``` r
de <- data.frame(
    logFC = c(2.1, -1.4, 0.2), logCPM = c(5, 7, 3),
    PValue = c(1e-8, 1e-4, 0.6), FDR = c(3e-8, 2e-4, 0.6),
    row.names = c("Cd3e", "Ms4a1", "Actb")
)
standardize_de(de)
#> Detected edgeR results (logFC, PValue, FDR)
#>   feature logfc pvalue  padj mean_expr
#> 1    Cd3e   2.1  1e-08 3e-08         5
#> 2   Ms4a1  -1.4  1e-04 2e-04         7
#> 3    Actb   0.2  6e-01 6e-01         3
standardize_de(list("KO vs WT" = de, "Het vs WT" = de))
#> Detected edgeR results in 'KO vs WT' (logFC, PValue, FDR)
#> Detected edgeR results in 'Het vs WT' (logFC, PValue, FDR)
#>   feature logfc pvalue  padj mean_expr comparison
#> 1    Cd3e   2.1  1e-08 3e-08         5   KO vs WT
#> 2   Ms4a1  -1.4  1e-04 2e-04         7   KO vs WT
#> 3    Actb   0.2  6e-01 6e-01         3   KO vs WT
#> 4    Cd3e   2.1  1e-08 3e-08         5  Het vs WT
#> 5   Ms4a1  -1.4  1e-04 2e-04         7  Het vs WT
#> 6    Actb   0.2  6e-01 6e-01         3  Het vs WT
```
