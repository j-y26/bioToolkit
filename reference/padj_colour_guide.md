# Breaks, labels, title and direction of an adjusted p-value colour bar

Breaks, labels, title and direction of an adjusted p-value colour bar

## Usage

``` r
padj_colour_guide(limits, padj.scale)
```

## Arguments

- limits:

  Limits of the colour scale, in the values of
  [`padj_colour_value()`](https://j-y26.github.io/bioToolkit/reference/padj_colour_value.md)

- padj.scale:

  "log", "-log10" or "linear"

## Value

A list: \`breaks\` and \`labels\` (waiver() for "-log10"), \`title\`,
and \`reverse\` (TRUE when the colour bar is flipped so that the most
significant values are at its top, as on the -log10 scales)
