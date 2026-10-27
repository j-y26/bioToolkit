# Write the negative numbers of a plot's scales with a true minus sign

ggplot2 labels negative breaks with the ASCII hyphen-minus ("-2"), which
is shorter than a minus and sits lower. The labels of every continuous
and binned scale (axes, colour bars, size legends, ...) are wrapped so
that a leading "-" before a number becomes U+2212. Scales added by
ggplot2 by default are found by building the plot once; the plot is then
rebuilt with them, so its layout is measured with the minus. Labels that
are not plain text (plotmath expressions) and discrete scales are left
as they are.

## Usage

``` r
minus_signs(plot)
```

## Arguments

- plot:

  A ggplot or patchwork object; anything else (spacers, grobs) is
  returned unchanged

## Value

The plot, with the labels of its continuous scales wrapped
