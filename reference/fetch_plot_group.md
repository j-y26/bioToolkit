# A grouping variable of the cells, from a Seurat object or a data frame

As
[`fetch_group()`](https://j-y26.github.io/bioToolkit/reference/fetch_group.md):
levels come from all cells (rows), not only from \`cells\`, so that
colours match a plot of all cells.

## Usage

``` r
fetch_plot_group(object, var, cells)
```

## Arguments

- object:

  A Seurat object, or a data frame with one row per cell

- var:

  Column name. NULL gives the active identities of a Seurat object, or
  one group ("All cells") for a data frame

- cells:

  Cell names or row indices

## Value

A factor along \`cells\`
