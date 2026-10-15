# A ggplot of fixed width holding a legend block at its top left

Only the width is fixed: placed in a patchwork over several rows of
fixed-size panels, the figure takes their height.

## Usage

``` r
legend_block_figure(block, pad = 0.15)
```

## Arguments

- block:

  A legend block
  ([`legend_block_grob()`](https://j-y26.github.io/bioToolkit/reference/legend_block_grob.md))

- pad:

  Room left of the legends, in inches

## Value

A ggplot object, with the block's measured height (inches, PDF metrics
plus 5
