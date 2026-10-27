# Least-squares trend line of a quadrant plot, with its band

Fits \`y ~ x\` and predicts it on 100 points across the range of \`x\`.
The band is clamped to the y range of the panel (a band beyond the scale
limits would be dropped by ggplot2); the line is drawn only where it
lies inside it.

## Usage

``` r
quadrant_trend(x, y, interval, ylim)
```

## Arguments

- x, y:

  Fold changes on x and y

- interval:

  "confidence", "prediction" or "none" (a band is still computed, as
  confidence, but not drawn)

- ylim:

  The y limits of the panel

## Value

A data frame with \`x\`, \`fit\`, \`lwr\`, \`upr\` and \`in_panel\`
