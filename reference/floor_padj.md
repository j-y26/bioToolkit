# Adjusted p-values of 0 replaced by the smallest non-zero one

Adjusted p-values of 0 replaced by the smallest non-zero one

## Usage

``` r
floor_padj(padj)
```

## Arguments

- padj:

  Adjusted p-values

## Value

\`padj\`, without zeros

## Details

Some tools report adjusted p-values of 0 (below machine precision),
whose -log10 is infinite. They are drawn as the most significant
non-zero value of the results (or the smallest positive double, when all
are 0).
