# Ridges of one feature

Ridges of one feature

## Usage

``` r
ridge_panel(d, group_cols, log, adjust, overlap, annotate)
```

## Arguments

- d:

  Long table of one feature
  ([`distribution_data()`](https://j-y26.github.io/bioToolkit/reference/distribution_data.md))

- group_cols:

  Colours named by group

- log:

  Whether the value axis is log10 (the bandwidth is then chosen on that
  scale)

- adjust:

  Bandwidth multiplier

- overlap:

  Height of the tallest ridge, in rows

- annotate:

  Whether to print the % of cells above 0 per group

## Value

A ggplot object (no scales for the value axis yet)
