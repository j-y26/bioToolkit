# Build the dimension arrow grob

Returns an empty grob tree of class "bt_dimension_arrows" that holds the
arrow geometry; the arrows and labels themselves are built when it is
drawn (`makeContent.bt_dimension_arrows()`), scaled with the size of the
viewport (i.e. panel) it is drawn in, and measured from its bottom left
corner.

## Usage

``` r
dimension_arrow_grob(
  reduction = "umap",
  dims = NULL,
  arrow_scale = 1,
  arrow_length = 2,
  x_coord_adjust = 0.25,
  y_coord_adjust = x_coord_adjust,
  arrow_size = 0.25,
  arrow_width = 0.2,
  arrow_line_width = 3,
  arrow_color = "black",
  text_font_size = 12,
  text_dst = 0.3
)
```

## Arguments

- reduction, dims:

  As in
  [`add_dimension_arrows`](https://j-y26.github.io/bioToolkit/reference/add_dimension_arrows.md)

- arrow_scale, arrow_length, x_coord_adjust, y_coord_adjust, arrow_size,
  arrow_width, arrow_line_width, arrow_color, text_font_size, text_dst:

  The arrow geometry, as the dot.case arguments of
  [`add_dimension_arrows`](https://j-y26.github.io/bioToolkit/reference/add_dimension_arrows.md)
  (e.g. \`arrow_length\` is \`arrow.length\`): sizes on a 3.5-inch panel

## Value

A grob tree
