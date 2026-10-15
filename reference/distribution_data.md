# Values of features per cell, with their groups, in long form

Values of features per cell, with their groups, in long form

## Usage

``` r
distribution_data(
  object,
  features,
  group.by,
  split.by,
  cells,
  assay,
  slot,
  zeros,
  log
)
```

## Arguments

- object:

  A Seurat object, or a data frame with one row per cell

- features, group.by, split.by, cells, assay, slot, zeros, log:

  As in
  [`VlnPlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/VlnPlot_bioToolkit.md)

## Value

A list with \`long\` (one row per cell and feature: \`cell\`, \`group\`,
\`split\`, \`feature\`, \`value\`, \`drawn\`), \`features\` (those
found), \`is_gene\`, \`assay\` (or NULL), \`cells\` (a data frame of the
cells with their \`unit\`) and \`group_levels\` (a factor with all
levels of \`group.by\`, for the colours)
