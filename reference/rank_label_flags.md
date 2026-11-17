# Which features of a rank plot are named

Which features of a rank plot are named

## Usage

``` r
rank_label_flags(df, label, label.n, label.side, label.features)
```

## Arguments

- df:

  Table from
  [`rank_positions()`](https://j-y26.github.io/bioToolkit/reference/rank_positions.md)

- label:

  "top", "selected" or "none"

- label.n:

  Number of features per labelled end, in each split

- label.side:

  "first", "last" or "both"

- label.features:

  Features to name, or NULL

## Value

A logical vector along the rows of \`df\`
