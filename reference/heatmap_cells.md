# Cells of a heatmap with their group, block and annotations

Cells of a heatmap with their group, block and annotations

## Usage

``` r
heatmap_cells(object, group.by, split.by, annotations, cells)
```

## Arguments

- object:

  A Seurat object

- group.by, split.by, annotations, cells:

  As in
  [`DoHeatmap_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/DoHeatmap_bioToolkit.md)

## Value

A list with \`cells\` (a data frame: \`cell\`, \`group\`, \`split\` or
NA, \`block\`, and one column per annotation) and \`group_levels\` (a
factor with all levels of \`group.by\`, for the colours)
