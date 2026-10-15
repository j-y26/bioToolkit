# Arrows naming the directions of change, with counts, for a volcano plot

Two arrows at height \`y\`, pointing left (down-regulated) and right
(up-regulated) from the centre of the x range, each with its label above
and its count below.

## Usage

``` r
volcano_arrow_layers(
  xlim,
  y,
  labels,
  counts,
  cols,
  text.size,
  count.size,
  linewidth,
  head.length
)
```

## Arguments

- xlim:

  x range the arrows are laid out in (the axis limits)

- y:

  Height of the arrows (data units)

- labels:

  Labels of the down and up arrows

- counts:

  Counts of down- and up-regulated features

- cols:

  Colours of the down and up arrows

- text.size:

  Size of the labels (mm)

- count.size:

  Size of the counts (mm)

- linewidth:

  Line width of the arrows (ggplot2 units)

- head.length:

  Length of the arrow heads in points

## Value

A list of ggplot2 layers
