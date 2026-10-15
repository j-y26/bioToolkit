# Legend labels of significant and other results

Both cutoffs are spelled out, since a result can have a small adjusted
p-value and still not be significant for its fold change.

## Usage

``` r
significance_labels(padj.cutoff, logfc.cutoff = 0)
```

## Arguments

- padj.cutoff:

  Adjusted p-value cutoff

- logfc.cutoff:

  Absolute log2 fold change cutoff (0: none)

## Value

Two labels: significant, not significant
