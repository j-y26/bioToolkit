# Flag the plotted cells that are highlighted

Flag the plotted cells that are highlighted

## Usage

``` r
highlight_flags(object, cells, cells.highlight)
```

## Arguments

- object:

  A Seurat object

- cells:

  The plotted cells

- cells.highlight:

  Cells to highlight as names or indices of the object's cells, or NULL
  to highlight all

## Value

A logical vector, one value per cell in \`cells\`
