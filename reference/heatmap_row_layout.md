# Positions and sizes of the rows of a heatmap

Positions and sizes of the rows of a heatmap

## Usage

``` r
heatmap_row_layout(ft, cell.height, panel.height)
```

## Arguments

- ft:

  Feature table from
  [`dot_feature_table()`](https://j-y26.github.io/bioToolkit/reference/dot_feature_table.md)

- cell.height:

  Height of one row (inches)

- panel.height:

  Height of all rows (inches), or NULL: rows of \`cell.height\`, up to 6
  in

## Value

A list with \`rows\` (\`ft\` with \`y\`, 1 at the bottom, and
\`block\`), \`heights\` (inches per block, top to bottom),
\`row.height\` and \`blocks\`
