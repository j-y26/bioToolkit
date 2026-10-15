# A metadata column for a set of cells, as a factor with all of its levels

Levels come from all cells of the object, not only from \`cells\`, so
that groups get the same colours in a plot of a subset of cells as in a
plot of all of them. Levels of a factor are kept; other values are
sorted.

## Usage

``` r
fetch_group(object, group.by, cells)
```

## Arguments

- object:

  A Seurat object

- group.by:

  Metadata column, or NULL for the active identities

- cells:

  Cell names

## Value

A factor along \`cells\`
