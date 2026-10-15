# Draw an arrow grob in one or all panels of a plot

\`annotation_custom()\` has no data of its own and is therefore drawn in
every panel. Passing \`data\` a single row that carries the facetting
variable restricts the grob to the matching panel, which is how the
arrows are kept to the bottom left panel of a split plot.

## Usage

``` r
apply_arrow_grob(plot, grob, data = NULL)
```

## Arguments

- plot:

  A ggplot object

- grob:

  The grob to draw, e.g. from
  [`dimension_arrow_grob`](https://j-y26.github.io/bioToolkit/reference/dimension_arrow_grob.md)

- data:

  A one-row data frame naming the panel to draw in, or NULL (default) to
  draw in every panel

## Value

a ggplot object
