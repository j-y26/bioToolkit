# Draw the panels sharing one colour scale as a single ggplot

Draw the panels sharing one colour scale as a single ggplot

## Usage

``` r
feature_ggplot(
  members,
  colour_scale,
  strip,
  strip.face,
  facet_ncol,
  facet_nrow,
  facet_dir,
  geom,
  order,
  center,
  pt.size,
  alpha,
  outline,
  outline.color,
  outline.size,
  label,
  label.size,
  label.color,
  legend.position,
  legend.direction,
  legend.length,
  legend.width,
  legend.text.size,
  legend.title.face,
  panel.width,
  panel.height,
  panel.spacing,
  strip.text.size,
  rd,
  raster,
  raster.threshold,
  raster.dpi,
  bg,
  halo.color,
  seed,
  arrow.scale = 1
)
```

## Arguments

- members:

  Panels to draw, each with its \`cells\`, \`marks\` and \`label\`

- colour_scale:

  Resolved colour scale, from
  [`resolve_colour_scale`](https://j-y26.github.io/bioToolkit/reference/resolve_colour_scale.md)

- strip:

  Whether to label the panels with a strip

- strip.face:

  Font face of the strip text

- facet_ncol, facet_nrow, facet_dir:

  Facet layout

- geom:

  How values are drawn: "point" (default) colours each cell; "hex"
  colours hexagonal tiles by a summary of the cells in them (see
  \`hex.fun\`); "density" colours each cell by values smoothed over
  neighbouring cells (see \`density.method\`)

- order:

  Whether to draw cells with higher values on top, default is TRUE; on a
  diverging scale, cells further from the centre are drawn on top. If
  FALSE, cells are drawn in random order

- center:

  Centre of a diverging scale, default is 0. The two sides of the scale
  extend equally far in colour from the centre, so that a side with a
  smaller range does not reach its darkest colour

- pt.size:

  Point size, default is 0.3

- alpha:

  Point (or tile) transparency, default is 1

- outline:

  Whether to draw an outline around the cell clouds, default is FALSE.
  Not used for \`geom = "hex"\`

- outline.color:

  Colour of the outline, default is "black"

- outline.size:

  Thickness of the outline (in point size units), default is 0.6

- label:

  Whether to label groups at their median position, default is FALSE.
  Groups are the active identities unless \`label.by\` is given

- label.size:

  Size of the group labels in mm, default is 3.5

- label.color:

  Colour of the group labels, default is "black"

- legend.position:

  Position of the legend, one of "right" (default), "left", "top",
  "bottom", "inside" or "none"

- legend.direction:

  Direction of the colour bar, "vertical" or "horizontal". If NULL
  (default), horizontal for top and bottom legends and vertical
  otherwise

- legend.length:

  Length of the colour bar in inches, default is 1.5

- legend.width:

  Thickness of the colour bar in inches, default is 0.15

- legend.text.size:

  Font size of the legend text and title in points. If NULL (default),
  the theme's base size (11)

- legend.title.face:

  Font face of the legend title. "auto" (default) sets it in bold, gene
  names in bold italics

- panel.width:

  Width of a single panel in inches, default is 4

- panel.height:

  Height of a single panel in inches, default is \`panel.width\`

- panel.spacing:

  Gap between panels in inches, default is 0.1

- strip.text.size:

  Font size of the panel labels in points. If NULL (default), the
  theme's large size (12)

- rd:

  Output of
  [`reduction_data`](https://j-y26.github.io/bioToolkit/reference/reduction_data.md)

- raster:

  Whether to rasterize the points while keeping text, legend and arrows
  as vectors: NULL (default) rasterizes when there are more than
  \`raster.threshold\` cells, or TRUE / FALSE. Hex tiles are always
  drawn as vectors

- raster.threshold:

  Number of cells above which \`raster = NULL\` rasterizes, default is
  10000

- raster.dpi:

  Resolution of the rasterized points, default is 1000. Match the
  \`dpi\` used when saving raster (e.g. png) output

- bg:

  Background colour, default is "transparent", so that the figure can be
  placed on any background (gaps inside outlines and label halos are
  then white). Use e.g. "white" for an opaque background

- halo.color:

  Colour of the gap inside outlines and of label halos

- seed:

  Random seed used for drawing order and label placement, default is 42.
  The global random number generator state is not changed

- arrow.scale:

  The \`arrow.scale\` of the dimension arrows, which sets the room left
  for their labels

## Value

a ggplot object
