# Hexagonal tiles for \`geom = "hex"\`

The hexagon grid spans the axis limits and is shaped to the panel, so
that hexagons are regular on the page and identical across panels.

## Usage

``` r
hex_marks(cells, v, cut, rd, bins, shape, fun, min_cells)
```

## Arguments

- cells:

  Cells of a panel, with x, y and highlight columns

- v:

  Values of the feature for these cells

- cut:

  The resolved \`na.cutoff\`

- rd:

  Output of
  [`reduction_data`](https://j-y26.github.io/bioToolkit/reference/reduction_data.md)

- bins:

  Number of hexagons across the width of a panel

- shape:

  Height to width ratio of a panel

- fun:

  "mean" or "pct"

- min_cells:

  Tiles with fewer plotted cells are not drawn

## Value

A data frame of hexagon vertices, six rows per tile, with tile, x, y,
value and coloured columns
