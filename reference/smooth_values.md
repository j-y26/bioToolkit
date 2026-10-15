# Kernel-smoothed values at the cells' own positions

The cells' values (and, for "smooth", their counts) are smoothed on a
grid
([`kde_grid()`](https://j-y26.github.io/bioToolkit/reference/kde_grid.md))
and read back at every cell by bilinear interpolation
([`grid_interpolate()`](https://j-y26.github.io/bioToolkit/reference/grid_interpolate.md)).

## Usage

``` r
smooth_values(x, y, w, xlim, ylim, bandwidth, method, n_grid = 256L)
```

## Arguments

- x, y:

  Coordinates of the cells

- w:

  Values of the cells

- xlim, ylim:

  Axis limits, which the grid spans

- bandwidth:

  Kernel standard deviations along x and y

- method:

  "smooth" for the kernel-weighted mean of \`w\`; "wkde" for the kernel
  density weighted by \`w\`, per cell

- n_grid:

  Number of grid points along each axis, default is 256

## Value

A numeric vector, one value per cell
