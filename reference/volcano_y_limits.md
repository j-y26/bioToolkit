# y-axis limits of a volcano plot

The panel has a fixed height, so a given space in inches above the
points can be converted to data units.

## Usage

``` r
volcano_y_limits(y.max, panel.height, headroom)
```

## Arguments

- y.max:

  Highest -log10 adjusted p-value to show

- panel.height:

  Panel height in inches

- headroom:

  Space in inches to keep free above the points

## Value

A list with the axis \`limits\`, the top of the points region
(\`points.top\`, the upper limit for labels) and the height of the
arrows (\`arrow.y\`, the middle of the headroom)
