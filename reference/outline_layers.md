# Outline layers drawn beneath all cells

A coloured layer then a background-coloured layer slightly larger than
the points, drawn beneath all cells, give the cell clouds a crisp
border.

## Usage

``` r
outline_layers(df, pt.size, outline.size, outline.color, halo.color, layer_fun)
```

## Arguments

- df:

  Plotted cells

- pt.size:

  Point size

- outline.size:

  Thickness of the outline

- outline.color:

  A single colour, or one colour per row of \`df\`

- halo.color:

  Colour of the gap between the outline and the cells

- layer_fun:

  Function wrapping point layers, from
  [`points_layer_fun`](https://j-y26.github.io/bioToolkit/reference/points_layer_fun.md)

## Value

A list of two layers
