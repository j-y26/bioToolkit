# Map features of a DE table (e.g. peaks) to genes, one per gene

Map features of a DE table (e.g. peaks) to genes, one per gene

## Usage

``` r
map_features(d, map, summary)
```

## Arguments

- d:

  Standardized DE table with \`feature\`, \`logfc\`, \`padj\`

- map:

  A data frame (first column: features, second: genes) or genes named by
  feature

- summary:

  "min_padj" or "max_abs"

## Value

\`d\` with one row per gene: \`feature\` is the gene, \`feature_y\` the
feature chosen and \`n_mapped\` the number of features of the gene
