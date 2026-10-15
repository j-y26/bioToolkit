# Resolve the \`cols\` argument of a continuous colour scale into colours

Shared by every continuous scale of the package (feature plots, dot
plots, heatmaps, density differences), so that they all accept the same
values.

## Usage

``` r
resolve_continuous_palette(cols, diverging, default = "plasma_r")
```

## Arguments

- cols:

  NULL for the default (\`default\`, or "blue-magenta" for a diverging
  scale); a sequential or diverging palette name of the package, or an
  RColorBrewer palette name; or two or more colours from low to high,
  e.g. \`c(low, high)\` or \`c(low, mid, high)\`

- diverging:

  Whether the scale is diverging, which sets the default

- default:

  Default sequential palette, "plasma_r" unless the plot sets another

## Value

A vector of colours, from low to high
