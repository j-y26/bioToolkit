# Normalized kernel density of one condition's cells on a grid

Normalized kernel density of one condition's cells on a grid

## Usage

``` r
condition_density(d, rd, bandwidth, n.grid, per_sample)
```

## Arguments

- d:

  Cells of the condition (\`x\`, \`y\` and, per sample, \`sample\`)

- rd:

  Result of
  [`reduction_data()`](https://j-y26.github.io/bioToolkit/reference/reduction_data.md)
  (axis limits)

- bandwidth:

  Kernel standard deviations

- n.grid:

  Grid size

- per_sample:

  Whether to average one density per sample

## Value

A grid summing to 1, with attributes \`n_cells\` and \`n_samples\`
