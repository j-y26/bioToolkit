# Add dimension arrows to one or all panels

Add dimension arrows to one or all panels

## Usage

``` r
add_reduction_arrows(p, reduction, labels, arrow_args, panel = NULL)
```

## Arguments

- p:

  A ggplot object

- reduction:

  Name of the reduction

- labels:

  The two arrow labels

- arrow_args:

  Further arguments to
  [`dimension_arrow_grob`](https://j-y26.github.io/bioToolkit/reference/dimension_arrow_grob.md)

- panel:

  Levels of the facetting variable \`panel\` to draw the arrows in, or
  NULL to draw them in every panel

## Value

a ggplot object
