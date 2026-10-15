# Embedding, axis limits and split panels for the plotted cells

Embedding, axis limits and split panels for the plotted cells

## Usage

``` r
reduction_data(object, reduction, dims, cells, split.by, full.extent)
```

## Arguments

- object:

  A Seurat object

- reduction:

  Name of the reduction

- dims:

  Two dimension indices

- cells:

  Cells to plot, or NULL for all

- split.by:

  Metadata column to split by, or NULL

- full.extent:

  Whether the axis limits span all cells of the reduction rather than
  only the plotted ones

## Value

A list with \`df\` (x, y and, when split, panel per plotted cell, with
cell names as row names), the axis limits \`xlim\` and \`ylim\`,
\`split_plot\` and the number of panels \`n_panels\`
