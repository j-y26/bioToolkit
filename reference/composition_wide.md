# One row per group with its values in the two conditions

One row per group with its values in the two conditions

## Usage

``` r
composition_wide(tab, x, y, value)
```

## Arguments

- tab:

  Result of
  [`summarize_group_meta()`](https://j-y26.github.io/bioToolkit/reference/summarize_group_meta.md)

- x, y:

  Conditions on the axes

- value:

  "percent" or "count"

## Value

A data frame: \`facet\` (if any), \`group\`, \`x_value\`, \`y_value\`,
and \`x_count\`, \`x_total\`, \`y_count\`, \`y_total\` (pooled), or
\`x_sd\`, \`x_sem\`, \`y_sd\`, \`y_sem\`, \`x_n_samples\`,
\`y_n_samples\` (per sample)
