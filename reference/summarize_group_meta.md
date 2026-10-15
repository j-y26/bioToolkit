# Count and proportion of cells per group and condition

Counts the cells of each group (e.g. cell type) in each condition, and
the proportion they make up of that condition's cells, optionally within
panels (\`facet.by\`) and per sample (\`sample.by\`). The data behind
[`plot_composition_scatter`](https://j-y26.github.io/bioToolkit/reference/plot_composition_scatter.md).

## Usage

``` r
summarize_group_meta(
  object,
  group.by,
  split.by,
  facet.by = NULL,
  sample.by = NULL,
  cells = NULL
)
```

## Arguments

- object:

  A Seurat object, or a data frame with one row per cell (e.g. its
  metadata)

- group.by:

  Column of the groups counted (e.g. cell type)

- split.by:

  Column of the conditions (e.g. genotype)

- facet.by:

  Column within whose levels proportions are computed separately (e.g.
  stage), or NULL (default)

- sample.by:

  Column of the samples, to count each sample separately, or NULL
  (default) to pool cells

- cells:

  Cells (rows) to use, as names or indices. If NULL (default), all

## Value

A data frame with one row per group and condition (and facet and
sample): \`facet\`, \`sample\` (when given), \`condition\`, \`group\`,
\`count\`, \`total\` (cells of the condition) and \`prop\`. Groups with
no cells in a condition are included with a count of 0. Cells without a
value in any of the columns are left out, with a message

## See also

[`plot_composition_scatter`](https://j-y26.github.io/bioToolkit/reference/plot_composition_scatter.md)

## Examples

``` r
meta <- data.frame(
    cell_type = c("A", "A", "B", "B", "B", "C"),
    genotype = c("WT", "Mut", "WT", "Mut", "Mut", "WT")
)
summarize_group_meta(meta, group.by = "cell_type", split.by = "genotype")
#>   condition group count total      prop
#> 1       Mut     A     1     3 0.3333333
#> 2       Mut     B     2     3 0.6666667
#> 3       Mut     C     0     3 0.0000000
#> 4        WT     A     1     3 0.3333333
#> 5        WT     B     1     3 0.3333333
#> 6        WT     C     1     3 0.3333333
```
