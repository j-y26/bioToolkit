# Combine one plot per feature

Combine one plot per feature

## Usage

``` r
combine_feature_plots(plots, ncol, stack, flip, bg)
```

## Arguments

- plots:

  Named list of ggplot objects

- ncol:

  Columns of the grid when not stacked, or NULL

- stack:

  Whether to stack them as rows (columns when \`flip\`)

- flip:

  Whether groups are on the y-axis

- bg:

  Background colour

## Value

The plot itself for one feature, otherwise a patchwork
