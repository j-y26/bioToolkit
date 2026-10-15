# One mark per cell, coloured by smoothed values, for \`geom = "density"\`

Only highlighted cells with a value take part in the smoothing, and only
they are coloured.

## Usage

``` r
density_marks(cells, v, cut, rd, bandwidth, method)
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

- bandwidth:

  Kernel standard deviations along x and y

- method:

  "smooth" or "wkde"
