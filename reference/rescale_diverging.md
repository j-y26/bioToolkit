# Rescaler of a diverging colour scale

Maps \`center\` to the middle of the palette and stretches each side of
it to its own end of the limits, for the \`rescaler\` of a ggplot2
continuous scale.

## Usage

``` r
rescale_diverging(center = 0)
```

## Arguments

- center:

  Centre of the scale

## Value

A function \`(x, to, from)\`
