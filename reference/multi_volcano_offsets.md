# Horizontal offsets of points within their column

"sina" (user decision 2026-10-15, after trying a shape per side of 0
built from all genes): the significant genes of each direction get their
own shape, as wide as their density, so their distribution shows however
many non-significant genes there are; non-significant genes are spread
uniformly over the column, a grey block centred on 0 rather than another
shape (built from all genes, a shape is dominated by the thousands of
non-significant genes and the significant ones become a thin line).
"pvalue" places genes by -log10 padj relative to the most significant
gene on their side of 0, significant or not, so that grey genes sit at
their own significance (the layout suited to bulk data, with many more
genes).

## Usage

``` r
multi_volcano_offsets(df, layout, half)
```

## Arguments

- df:

  Genes with \`comparison\`, \`direction\`, \`logfc\` and
  \`neg_log10_padj\`

- layout:

  "sina", "pvalue" or "jitter"

- half:

  Half the width of a column, in column units

## Value

Offsets in \[-half, half\], along the rows of \`df\`
