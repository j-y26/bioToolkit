# Smooth outline around the cloud of a set of cells

The cells are binned onto a grid and smoothed with a Gaussian kernel of
sd \`smooth\`. The region where this density is at least its value at
the \`1 - coverage\` quantile of the cells is the cloud: the sparsest
cells, such as strays far from the rest, fall outside it. Holes in the
cloud are filled. The cloud is then blurred and cut at the level that,
along a straight edge, lies exactly \`expand\` beyond it, which both
pushes the outline away from the cells and rounds its corners. Anything
enclosed by the outline is filled again so that only outer outlines are
returned.

## Usage

``` r
cloud_outline(u, v, expand, smooth, coverage, resolution = 0.01)
```

## Arguments

- u, v:

  Coordinates of the cells, in inches

- expand:

  Distance between the cloud and its outline, in inches

- smooth:

  Smoothing of the cloud, in inches. Larger values give rounder outlines
  that bridge larger gaps between cells

- coverage:

  Fraction of the cells, from the densest regions, taken to form the
  cloud

- resolution:

  Spacing of the grid, in inches

## Value

A data frame with u, v and ring columns, one closed ring per separate
part of the cloud
