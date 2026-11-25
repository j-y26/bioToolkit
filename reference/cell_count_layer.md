# Number of cells written in a corner of each panel

Number of cells written in a corner of each panel

## Usage

``` r
cell_count_layer(
  df,
  by_panel,
  xlim,
  ylim,
  position = "bottom.right",
  size,
  inset = 0.02
)
```

## Arguments

- df:

  Plotted cells, with a \`panel\` column when \`by_panel\`

- by_panel:

  Whether to count cells within each panel

- xlim, ylim:

  Axis limits of the panels (data units)

- position:

  Corner of the panel: "bottom.right", "bottom.left", "top.right" or
  "top.left"

- size:

  Text size in ggplot2 size units (mm)

- inset:

  Distance from the panel edges, as a fraction of the axis ranges

## Value

A text layer with one label per panel
