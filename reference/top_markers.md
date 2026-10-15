# Top marker genes per group

Picks the top marker genes of each group from
[`Seurat::FindAllMarkers()`](https://satijalab.org/seurat/reference/FindAllMarkers.html)
output (or any DE table with a group column), as a named list for the
\`features\` of
[`DotPlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/DotPlot_bioToolkit.md).

## Usage

``` r
top_markers(
  markers,
  n = 5,
  group.col = "cluster",
  rank.by = c("logfc", "padj"),
  padj.cutoff = 0.05,
  unique = TRUE,
  ...
)
```

## Arguments

- markers:

  A data frame of markers, e.g. from
  [`Seurat::FindAllMarkers()`](https://satijalab.org/seurat/reference/FindAllMarkers.html);
  read with
  [`standardize_de`](https://j-y26.github.io/bioToolkit/reference/standardize_de.md)

- n:

  Number of genes per group, default is 5

- group.col:

  Column naming the group of each row, default is "cluster" (as
  FindAllMarkers)

- rank.by:

  How genes are ranked within a group: "logfc" (default; largest log2
  fold change, ties broken by the adjusted p-value) or "padj" (smallest
  adjusted p-value, ties broken by the largest log2 fold change).
  Ranking by p-value favours broadly expressed genes in large groups;
  ranking by fold change favours specific genes

- padj.cutoff:

  Only genes with an adjusted p-value below this, and a positive log2
  fold change, are markers. Default is 0.05

- unique:

  Whether each gene is a marker of one group only, default is TRUE: a
  gene that is a marker of several groups is kept for the group where
  its log2 fold change is largest

- ...:

  Further arguments to
  [`standardize_de`](https://j-y26.github.io/bioToolkit/reference/standardize_de.md),
  e.g. \`feature.col\`

## Value

A named list of character vectors, one per group with markers, in the
order of the groups

## See also

[`DotPlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/DotPlot_bioToolkit.md)

## Examples

``` r
markers <- data.frame(
    gene = c("Cd3e", "Cd3d", "Ms4a1", "Cd79a", "Cd3e"),
    cluster = factor(c("T", "T", "B", "B", "B")),
    avg_log2FC = c(3, 2.5, 4, 3.5, 0.5),
    p_val = 1e-10, p_val_adj = 1e-8
)
top_markers(markers, n = 2)
#> $B
#> [1] "Ms4a1" "Cd79a"
#> 
#> $T
#> [1] "Cd3e" "Cd3d"
#> 
```
