# Values of a heatmap: z-scored cells (downsampled) or group means

Values of a heatmap: z-scored cells (downsampled) or group means

## Usage

``` r
heatmap_values(
  expr,
  hd,
  columns,
  average,
  scale,
  is_gene,
  log.data,
  max.cells,
  seed
)
```

## Arguments

- expr:

  Matrix of cells x features, rows in the order of \`hd\$cells\`

- hd:

  Result of
  [`heatmap_cells()`](https://j-y26.github.io/bioToolkit/reference/heatmap_cells.md)

- columns:

  "cells" or "groups"

- average:

  "linear" or "log"

- scale:

  Whether to z-score each feature

- is_gene:

  Logical, named by feature

- log.data:

  Whether genes hold log1p values (the "data" layer)

- max.cells:

  Largest number of cells per group and block, or NULL

- seed:

  Random seed

## Value

A list with \`columns\` (a data frame, one row per column: \`column\`,
\`block\`, \`group\`, \`split\`, annotations and \`n_cells\`), \`value\`
and \`fill\` (matrices of columns x features: values as read or
averaged, and as coloured)
