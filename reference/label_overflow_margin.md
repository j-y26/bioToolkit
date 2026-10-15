# Left margin for a first group label that runs past the plot's left edge

Group labels of vertical violins are turned 45 degrees and end under
their violin, so a long first label reaches left of the panel, past the
axis and out of the saved figure. The margin is estimated from the
number of characters (about 0.55 em per character) less the room left of
the violin (half a slot and about 0.3 in for the axis).

## Usage

``` r
label_overflow_margin(first, cell.width, sz, margins)
```

## Arguments

- first:

  First group label

- cell.width:

  Width of one violin's slot

- sz:

  Sizes from
  [`theme_sizes()`](https://j-y26.github.io/bioToolkit/reference/theme_sizes.md)

- margins:

  Plot margins (top, right, bottom, left) in points

## Value

The margins, with the left one widened if needed
