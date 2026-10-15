# Grid layout of a set of panels

Grid layout of a set of panels

## Usage

``` r
panel_layout(n, ncol, by.row, default_ncol)
```

## Arguments

- n:

  Number of panels

- ncol:

  Requested number of columns, or NULL for \`default_ncol\`

- by.row:

  Whether panels fill the grid row by row (TRUE) or column by column
  (FALSE)

- default_ncol:

  Number of columns used when \`ncol\` is NULL

## Value

A list with \`ncol\`, \`nrow\`, the facet direction \`dir\` and the
index of the bottom left panel, \`bottom_left\`
