# Smoothing bandwidth for \`geom = "density"\`

Scott's rule in two dimensions, with a robust estimate of the spread so
that the many separate clusters of an embedding do not inflate it.

## Usage

``` r
density_bandwidth(x, y, adjust)
```

## Arguments

- x, y:

  Coordinates of all plotted cells

- adjust:

  Multiplier of the bandwidth

## Value

Kernel standard deviations along x and y
