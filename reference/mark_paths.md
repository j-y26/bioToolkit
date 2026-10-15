# Outlines around the cell clouds of each group, in data units

Cells are converted to inches on the panel, outlined by
[`cloud_outline`](https://j-y26.github.io/bioToolkit/reference/cloud_outline.md)
per group (and per panel, when split), and converted back, so that the
outlines keep the same distance from the cells and the same smoothness
on the page whatever the range of the embedding.

## Usage

``` r
mark_paths(
  df,
  by_panel,
  xlim,
  ylim,
  panel.width,
  panel.height,
  expand,
  smooth,
  coverage
)
```

## Arguments

- df:

  Cells with x, y, mark and, when \`by_panel\`, panel columns

- by_panel:

  Whether to outline groups within each panel

- xlim, ylim:

  Axis limits of the panel

- panel.width, panel.height:

  Size of a panel in inches

- expand:

  Distance between the cloud and its outline, in inches

- smooth:

  Smoothing of the cloud, in inches. Larger values give rounder outlines
  that bridge larger gaps between cells

- coverage:

  Fraction of the cells, from the densest regions, taken to form the
  cloud

## Value

A data frame with x, y, mark, ring (one per closed outline) and, when
\`by_panel\`, panel columns
