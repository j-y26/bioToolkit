# Violins, points and boxes of one feature

Violins, points and boxes of one feature

## Usage

``` r
violin_panel(
  d,
  pts,
  group_cols,
  flip,
  box,
  box.width,
  pt.size,
  alpha,
  adjust,
  layer_fun,
  seed,
  annotate,
  top = NA_real_,
  split_side = if (flip) "right" else "top"
)
```

## Arguments

- d:

  Long table of one feature
  ([`distribution_data()`](https://j-y26.github.io/bioToolkit/reference/distribution_data.md))

- pts:

  Rows of \`d\` drawn as points

- group_cols:

  Colours named by group

- flip:

  Whether violins are horizontal

- box, box.width, pt.size, alpha, adjust, seed:

  As in
  [`VlnPlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/VlnPlot_bioToolkit.md)

- layer_fun:

  Function rasterizing the point layer, or not

- annotate:

  Whether to print the % of cells above 0 per group

- top:

  Top of the value axis on the scale it is drawn on (log10 for a log
  axis), where violins and boxes are cut; NA for none

- split_side:

  Side of the split.by block names: "top" / "bottom" (vertical violins)
  or "left" / "right" (horizontal)

## Value

A ggplot object (no scales for the value axis yet)
