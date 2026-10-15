# Which features of a volcano plot are labelled

Which features of a volcano plot are labelled

## Usage

``` r
volcano_label_flags(df, direction, label, label.n, label.features, label.rank)
```

## Arguments

- df:

  Standardized DE table with \`feature\`, \`logfc\`, \`padj\`

- direction:

  "up", "down" or "ns" per row of \`df\`

- label:

  "top", "selected" or "none"

- label.n:

  Number of top features per direction

- label.features:

  Features to label, or NULL

- label.rank:

  "padj", "logfc" or "product"

## Value

A logical vector along the rows of \`df\`
