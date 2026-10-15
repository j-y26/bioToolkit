# The plots, layout and annotation of a patchwork

A patchwork is the last plot added to it, carrying the earlier plots in
its \`patches\`; the last plot is split off and appended to them, unless
it is the empty placeholder that
[`patchwork::wrap_plots()`](https://patchwork.data-imaginist.com/reference/wrap_plots.html)
creates.

## Usage

``` r
patchwork_parts(plot)
```

## Arguments

- plot:

  A patchwork object

## Value

A list with \`plots\`, \`layout\` and \`annotation\`
