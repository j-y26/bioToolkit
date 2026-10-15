# Build a plot and take out its legend

Builds the plot on the current device. A hidden legend is drawn on the
right instead.

## Usage

``` r
legend_gtable(plot)
```

## Arguments

- plot:

  A ggplot

## Value

A list with the \`legend\` (a gtable), the full plot gtable (\`gt\`) and
the \`plot\` it was built from (with the legend shown)
