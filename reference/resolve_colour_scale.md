# Resolve the colour scale shared by a group of panels

Resolve the colour scale shared by a group of panels

## Usage

``` r
resolve_colour_scale(
  members,
  scale_features,
  min.cutoff,
  max.cutoff,
  diverging,
  center,
  cols,
  na.value,
  key
)
```

## Arguments

- members:

  Panels sharing the scale

- scale_features:

  Features of these panels

- min.cutoff, max.cutoff:

  Lists of cutoffs named by feature

- diverging:

  Whether to use a diverging scale centred on \`center\`: "auto"
  (default) does so when the coloured values span \`center\`, or TRUE /
  FALSE

- center:

  Centre of a diverging scale, default is 0. The two sides of the scale
  extend equally far in colour from the centre, so that a side with a
  smaller range does not reach its darkest colour

- cols:

  The colours of the scale. If NULL (default), "plasma_r" or, for a
  diverging scale, "blue-magenta". Otherwise two or more colours from
  low to high values, e.g. \`c("grey90", "darkred")\` (low, high) or,
  for a diverging scale, \`c("#2166AC", "white", "#B2182B")\` (low,
  centre, high); or the name of a sequential or diverging palette (see
  [`show_palettes`](https://j-y26.github.io/bioToolkit/reference/show_palettes.md)),
  such as "viridis", "magma" or "tol-burd", or of an RColorBrewer
  palette

- na.value:

  Colour of cells (or tiles) that are not coloured. If NULL (default),
  "grey88", or "grey70" for a diverging scale so that they stand out
  from its light centre

- key:

  Name of the scale, for messages

## Value

A list with the colour \`limits\`, \`colours\`, \`rescaler\`,
\`na.value\` and whether it is \`diverging\`
