# Resolve a cutoff given as a number, a quantile ("q##") or NA

Resolve a cutoff given as a number, a quantile ("q##") or NA

## Usage

``` r
cutoff_value(cutoff, v, default)
```

## Arguments

- cutoff:

  The cutoff

- v:

  Values the quantile or default is computed from

- default:

  Function giving the end when \`cutoff\` is NA (min or max)

## Value

A number, NA when \`v\` is empty
