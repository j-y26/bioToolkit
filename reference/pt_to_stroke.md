# Convert a line width in points to the \`stroke\` of ggplot2 points

ggplot2 draws the outline of a point \`stroke \* ggplot2::.stroke / 2\`
device line widths wide, each 1/96 inch, i.e. 0.75 pt.

## Usage

``` r
pt_to_stroke(pt)
```

## Arguments

- pt:

  Outline width in points

## Value

The \`stroke\` giving that width
