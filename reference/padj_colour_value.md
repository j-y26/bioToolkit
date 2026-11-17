# Colour values of adjusted p-values

Colour values of adjusted p-values

## Usage

``` r
padj_colour_value(padj, padj.scale)
```

## Arguments

- padj:

  Adjusted p-values

- padj.scale:

  "log", "-log10" or "linear"

## Value

-log10 padj (zeros floored, see
[`floor_padj()`](https://j-y26.github.io/bioToolkit/reference/floor_padj.md))
for "log" and "-log10", padj for "linear"
