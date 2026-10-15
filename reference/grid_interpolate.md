# Values of a grid at given positions, by bilinear interpolation

Values of a grid at given positions, by bilinear interpolation

## Usage

``` r
grid_interpolate(grid, x, y, xlim, ylim)
```

## Arguments

- grid:

  A square matrix from
  [`kde_grid()`](https://j-y26.github.io/bioToolkit/reference/kde_grid.md)

- x, y:

  Positions

- xlim, ylim:

  Axis limits the grid spans

## Value

A numeric vector along \`x\`
