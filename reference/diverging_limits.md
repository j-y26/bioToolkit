# Limits of a diverging colour scale, each end trimmed to the data

Limits of a diverging colour scale, each end trimmed to the data

## Usage

``` r
diverging_limits(v, low = NA, high = NA, center = 0)
```

## Arguments

- v:

  Values coloured

- low, high:

  Cutoffs of the two ends (NA for none): the end is the cutoff, or the
  data's extreme on that side when it is closer to the centre

- center:

  Centre of the scale

## Value

Two numbers, the low end at most \`center\` and the high end at least
\`center\`
