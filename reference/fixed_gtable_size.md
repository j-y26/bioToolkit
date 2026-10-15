# Absolute size of a gtable, including legends that overflow their cell

A legend taller than a fixed-height panel (right/left) or wider than a
fixed-width panel (top/bottom) spills outside its cell and is not
reflected in the gtable's widths/heights. The legend is centred on its
cell and the gtable is centred on the device, so the device is enlarged
just enough for the legend to keep \`margin\` inches from the device
edge
([`overflow_fit()`](https://j-y26.github.io/bioToolkit/reference/overflow_fit.md));
a legend that spills into the strips or titles above or below its cell,
but stays within the gtable, needs no more room. (Before bioToolkit
0.7.0, single plots were enlarged to the legend size plus twice the
larger of the spaces before and after its cell, which added up to
several inches of blank space to grid plots with tall strips.)

## Usage

``` r
fixed_gtable_size(gt, patchwork = FALSE, margin = 0.1)
```

## Arguments

- gt:

  A gtable with absolute widths and heights

- patchwork:

  Whether \`gt\` is the gtable of a patchwork

- margin:

  Smallest gap between a legend and the device edge, in inches

## Value

Named numeric vector with width and height in inches
