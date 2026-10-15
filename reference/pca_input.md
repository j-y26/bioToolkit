# PCA of a matrix or a \`prcomp\` result

PCA of a matrix or a \`prcomp\` result

## Usage

``` r
pca_input(x, n.top, center, scale)
```

## Arguments

- x:

  A \`prcomp\` result, or a features x samples matrix / data frame

- n.top:

  Most variable features to keep (NULL or Inf for all)

- center, scale:

  Passed to \`prcomp()\`

## Value

A list with \`scores\` (samples x PCs, columns "PC1", ...) and
\`variance\` (share of the total variance per PC)
