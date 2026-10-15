# Turn block (split) labels upright when they are wider than their block

Widths are estimated from the number of characters (bold text, about 0.6
em per character), as in the dot grid.

## Usage

``` r
block_label_theme(d, flip, cell.width, sz, side = "top", bar.width = 3.5)
```

## Arguments

- d:

  Long table of one feature

- flip:

  Whether blocks are rows (labels to the left, never turned)

- cell.width:

  Width of one violin's slot

- sz:

  Sizes from
  [`theme_sizes()`](https://j-y26.github.io/bioToolkit/reference/theme_sizes.md)

- side:

  Side of the block names ("top" or "bottom" when not flipped)

- bar.width:

  Thickness of the block bars in points

## Value

A theme
