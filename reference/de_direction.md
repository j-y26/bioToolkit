# Direction of change of DE features

A feature is up-regulated when its adjusted p-value is below
\`padj.cutoff\` and its log2 fold change above \`logfc.cutoff\`,
down-regulated when below \`-logfc.cutoff\`, and not significant
otherwise (also when the adjusted p-value is NA).

## Usage

``` r
de_direction(logfc, padj, padj.cutoff, logfc.cutoff)
```

## Arguments

- logfc:

  Log2 fold changes

- padj:

  Adjusted p-values

- padj.cutoff:

  Adjusted p-value cutoff

- logfc.cutoff:

  Absolute log2 fold change cutoff

## Value

"up", "down" or "ns" per feature
