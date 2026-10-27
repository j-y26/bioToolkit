# Which features of a quadrant plot are labelled

Which features of a quadrant plot are labelled

## Usage

``` r
quadrant_label_flags(df, both, label, label.n, label.features)
```

## Arguments

- df:

  Joined table with \`feature\`, \`logfc_x\`, \`logfc_y\`

- both:

  Whether each row is in a coloured class significant in both
  comparisons (the candidates of \`label = "top"\`)

- label:

  "top", "selected" or "none"

- label.n:

  Features per quadrant

- label.features:

  Features to label, or NULL

## Value

A logical vector along the rows of \`df\`
