# Panels of a grid that get the dimension arrows

Panels are numbered in the order they fill the grid (by row or by
column). The edges are those of the panels actually drawn: in a grid
whose last row is incomplete, a panel with no panel below it is on the
"bottom" edge, and the last panel of a row is on the "right" edge.

## Usage

``` r
arrow_panel_index(n, ncol, nrow, by.row, which)
```

## Arguments

- n:

  Number of panels

- ncol, nrow:

  Size of the grid

- by.row:

  Whether panels fill the grid row by row

- which:

  One of \`ARROW_PANELS\`: "bottom.left" (the first panel of the bottom
  row), "all", "left" / "right" (first / last panel of each row), "top"
  / "bottom" (top / lowest panel of each column) or "none"

## Value

Indices of the panels, in 1..n
