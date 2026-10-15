# Give a patchwork the background of its plots

A patchwork draws its own background, taken from ggplot2's global theme
(white by default) unless set with
[`patchwork::plot_annotation()`](https://patchwork.data-imaginist.com/reference/plot_annotation.html).
When it is not set, the background of the first plot is used, so that a
patchwork of transparent plots stays transparent and one of white plots
stays white.

## Usage

``` r
patchwork_background(plot)
```

## Arguments

- plot:

  A patchwork object

## Value

The patchwork, with a background set
