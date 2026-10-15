# Layer of group annotation bars

Layer of group annotation bars

## Usage

``` r
group_bar_layer(
  data,
  side,
  width = 3.5,
  offset = grid::unit(GROUP_BAR_GAP, "pt"),
  half = 0.5
)
```

## Arguments

- data:

  One row per bar: the facet columns of the panel it is drawn on, \`x\`
  / \`xend\` (bars on top or bottom) or \`y\` / \`yend\` (bars left or
  right) giving the first and last category of the group, and \`fill\`
  (colours)

- side:

  Side of the panel: "top", "bottom", "left" or "right"

- width:

  Thickness of the bar in points

- offset:

  Distance from the panel edge to the bar, a grid unit

- half:

  Half the width of a category, in data units (0.5 on a discrete scale)

## Value

A ggplot2 layer
