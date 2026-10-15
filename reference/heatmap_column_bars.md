# Bars of the blocks of columns of a heatmap, above the top row

Bars of the blocks of columns of a heatmap, above the top row

## Usage

``` r
heatmap_column_bars(cl, cols, width, row_block = NULL)
```

## Arguments

- cl:

  Result of
  [`heatmap_column_layout()`](https://j-y26.github.io/bioToolkit/reference/heatmap_column_layout.md)

- cols:

  Bar colours named by block

- width:

  Thickness of the bars in points

- row_block:

  When the bars are drawn on the body, its blocks of rows (the bars go
  on the first); NULL for a header row

## Value

A ggplot2 layer
