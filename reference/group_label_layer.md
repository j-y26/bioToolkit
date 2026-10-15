# Labels at the median position of each group

Labels at the median position of each group

## Usage

``` r
group_label_layer(df, by_panel, label.size, label.color, halo.color, seed)
```

## Arguments

- df:

  Cells with x, y, group and, when \`by_panel\`, panel columns

- by_panel:

  Whether to label groups within each panel

- label.size:

  Size of the labels

- label.color:

  Colour of the labels

- halo.color:

  Colour of the halo around the labels

- seed:

  Seed for label placement

## Value

A ggrepel text layer
