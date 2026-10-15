# Fetch one variable for a fixed set of cells

\`Seurat::FetchData()\` silently drops cells whose value is \`NA\`,
which would misalign the result with the embedding it is meant to
annotate. The values are therefore indexed back onto \`cells\`,
restoring the dropped ones as \`NA\` so that they can be shown as their
own group or panel.

## Usage

``` r
fetch_meta(object, var, cells)
```

## Arguments

- object:

  A Seurat object

- var:

  Name of the variable to fetch

- cells:

  Cell names to fetch, in the order they are wanted

## Value

A vector with one value per cell in \`cells\`
