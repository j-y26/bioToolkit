# Long table of a dot plot

Long table of a dot plot

## Usage

``` r
dot_long_table(ft, units, st, fill_mat, scale)
```

## Arguments

- ft:

  Gene table from
  [`dot_feature_table()`](https://j-y26.github.io/bioToolkit/reference/dot_feature_table.md)

- units:

  Unit table from
  [`dot_units()`](https://j-y26.github.io/bioToolkit/reference/dot_units.md)

- st:

  Statistics from
  [`group_expression_stats()`](https://j-y26.github.io/bioToolkit/reference/group_expression_stats.md)

- fill_mat:

  Values the dots are coloured by, units x genes

- scale:

  Whether \`fill_mat\` holds z-scores

## Value

A data frame, one row per unit and row of \`ft\`
