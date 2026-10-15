# Draw a grid of dots

Draw a grid of dots

## Usage

``` r
dot_grid_plot(
  df,
  x.labels,
  y.labels,
  x.italic = FALSE,
  y.italic = FALSE,
  size.limits,
  size.breaks,
  size.labels = waiver(),
  size.title,
  fill.colours,
  fill.limits,
  fill.center = NULL,
  fill.title,
  na.value = "grey85",
  outline.labels = NULL,
  outline.title = NULL,
  strip.axis = NULL,
  strip.style = "strip",
  cell.width = 0.18,
  cell.height = 0.18,
  x.group.side = "opposite",
  y.group.side = "opposite",
  x.group.cols = NULL,
  y.group.cols = NULL,
  bar.width = 3.5,
  legend.position = "right",
  bg = "transparent"
)
```

## Arguments

- df:

  One row per dot, with columns \`x\` and \`y\` (factors; their levels
  give the order of the columns from left to right and of the rows from
  top to bottom, and must be unique across groups), \`size\` (NA: no
  dot, but the row and column are kept), \`fill\`, and optionally
  \`x_group\` / \`y_group\` (factors: column / row groups), \`outline\`
  (logical: significant dots, drawn opaque with a dark ring; the others
  are faded) and \`strip\` (the colour of the row or column, for the
  group markers)

- x.labels, y.labels:

  Display labels named by level of \`x\` / \`y\`

- x.italic, y.italic:

  Whether each label is set in italics (genes): logical vectors named by
  level, or a single value

- size.limits:

  Range of \`size\` mapped to dot areas; values beyond are drawn at the
  nearest end

- size.breaks, size.labels:

  Breaks and labels of the size legend

- size.title:

  Title of the size legend

- fill.colours:

  Colours of the fill scale, from low to high

- fill.limits:

  Limits of the fill scale; values beyond get the colour of the nearest
  end. For a diverging scale they need not be symmetric (see
  [`diverging_limits()`](https://j-y26.github.io/bioToolkit/reference/diverging_limits.md))

- fill.center:

  Centre of a diverging fill scale, or NULL for a sequential one

- fill.title:

  Title of the fill legend

- na.value:

  Fill of dots whose \`fill\` is NA

- outline.labels:

  Legend labels of significant and other dots, when \`df\$outline\` is
  given

- outline.title:

  Title of the significance legend

- strip.axis:

  "x" when the group markers run along the columns (one per column,
  below the grid), "y" when they run along the rows (one per row, left
  of the grid), or NULL for none

- strip.style:

  "strip" (a narrow band of tiles) or "dot" (one coloured dot per row or
  column)

- cell.width, cell.height:

  Size of one column / row in inches

- x.group.side, y.group.side:

  Side of the column / row groups: "opposite" (default; opposite the
  tick labels: top / right) or "labels" (beyond the labels: bottom /
  left)

- x.group.cols, y.group.cols:

  Colours of the group bars, named by group; NULL for black

- bar.width:

  Thickness of the group bars in points

- legend.position:

  Position of the legends

- bg:

  Background colour

## Value

A ggplot object
