# Tidy table of a rank plot

Tidy table of a rank plot

## Usage

``` r
rank_input(data, rank.by, value, feature.col, split.by, color.by)
```

## Arguments

- data:

  Input data frame

- rank.by, value, feature.col, split.by, color.by:

  Columns, as in
  [`plot_rank`](https://j-y26.github.io/bioToolkit/reference/plot_rank.md)

## Value

A data frame with \`feature\`, \`split\` (factor, one level without
\`split.by\`), \`rank_value\`, \`value\` and \`color\` (NULL without
\`color.by\`), without the rows missing \`rank.by\` or \`value\`
