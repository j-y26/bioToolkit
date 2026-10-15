# Positions and sizes of the columns of a heatmap

Columns sit at integer positions within each block (facet). For cells,
groups within a block are set apart by empty columns about 0.015 in
wide.

## Usage

``` r
heatmap_column_layout(col_df, columns, panel.width, cell.width)
```

## Arguments

- col_df:

  Column table from
  [`heatmap_values()`](https://j-y26.github.io/bioToolkit/reference/heatmap_values.md)

- columns:

  "cells" or "groups"

- panel.width:

  Width of all columns of cells (inches)

- cell.width:

  Width of a column of group means (inches)

## Value

A list with \`columns\` (\`col_df\` with \`x\`), \`widths\` (inches per
block), \`width\` (one column, inches), \`limits\` (x range per block)
and \`blocks\` (block levels)
