# Scale expansion that makes categories tile a grid panel exactly

A discrete scale of n categories spans 1 to n; expanding it by half a
category on each side makes it span n units, so that each category gets
exactly the \`cell\` size given to
[`grid_panel_size`](https://j-y26.github.io/bioToolkit/reference/grid_panel_size.md).
(ggplot2's default of 0.6 would make cells slightly smaller than asked.)

## Usage

``` r
grid_expansion()
```

## Value

An expansion vector for the \`expand\` argument of a discrete scale
