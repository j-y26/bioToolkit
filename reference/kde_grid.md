# Gaussian kernel sum of weighted cells on a regular grid

Cells are binned onto a regular grid over the axis limits (nearest grid
point), and the binned weights are convolved with a Gaussian kernel,
applied along x and along y as two matrix products. With unit weights,
the result is a kernel density estimate up to a constant factor. Shared
by \`geom = "density"\` of
[`FeaturePlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/FeaturePlot_bioToolkit.md)
and by
[`DensityDiffPlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/DensityDiffPlot_bioToolkit.md).

## Usage

``` r
kde_grid(x, y, w, xlim, ylim, bandwidth, n_grid = 256L)
```

## Arguments

- x, y:

  Coordinates of the cells

- w:

  Weights of the cells

- xlim, ylim:

  Axis limits, which the grid spans

- bandwidth:

  Kernel standard deviations along x and y

- n_grid:

  Number of grid points along each axis

## Value

An \`n_grid\` x \`n_grid\` matrix: rows along x, columns along y
