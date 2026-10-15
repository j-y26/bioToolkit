# Colour bar breaks whose labels fit side by side along a bar

Breaks are chosen as usual, then every other one is dropped until their
labels, estimated from their number of characters, fit along the bar.

## Usage

``` r
fitted_breaks(length, text_size)
```

## Arguments

- length:

  Length of the bar in inches

- text_size:

  Font size of the labels

## Value

A function of the limits returning the breaks
