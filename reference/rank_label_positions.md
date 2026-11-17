# Positions of the labels of a rank plot

Each panel is laid out on its own, in inches, from its limits and the
panel size (see
[`rank_place_labels()`](https://j-y26.github.io/bioToolkit/reference/rank_place_labels.md)).
The layout is tried at each size of \`sizes\` in turn, until the labels
of every panel fit clear of the points and of each other, so that all
panels share one label size.

## Usage

``` r
rank_label_positions(
  df,
  lims,
  panel.width,
  panel.height,
  sizes,
  italic,
  decreasing
)
```

## Arguments

- df:

  Plotted table, with \`split\`, \`rank\`, \`value\`, \`labelled\`,
  \`.size\` and \`.label\`

- lims:

  Panel limits from
  [`rank_panel_limits()`](https://j-y26.github.io/bioToolkit/reference/rank_panel_limits.md)

- panel.width, panel.height:

  Panel size in inches

- sizes:

  Font sizes of the labels to try, in points, largest first

- italic:

  Whether the labels are italic

- decreasing:

  Whether rank 1 is the largest value (the curve falls)

## Value

The labelled rows of \`df\`, with \`label_x\`, \`label_y\` (centre of
the label) and \`seg_x\`, \`seg_y\` (end of the leader line), in data
units; attribute \`label.size\`, the size used
