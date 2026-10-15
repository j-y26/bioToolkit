# Long table of heatmap tiles, and the source data

Long table of heatmap tiles, and the source data

## Usage

``` r
heatmap_long(sv, ft, cl, rl)
```

## Arguments

- sv:

  Result of
  [`heatmap_values()`](https://j-y26.github.io/bioToolkit/reference/heatmap_values.md)

- ft:

  Feature table (rows of the heatmap)

- cl:

  Result of
  [`heatmap_column_layout()`](https://j-y26.github.io/bioToolkit/reference/heatmap_column_layout.md)

- rl:

  Result of
  [`heatmap_row_layout()`](https://j-y26.github.io/bioToolkit/reference/heatmap_row_layout.md)

## Value

A list with \`tiles\` (x, y, fill, block, row_block) and \`source\`
