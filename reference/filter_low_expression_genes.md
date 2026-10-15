# Filter out genes with low expression based on minimum group size

This function filters out genes with low expression based on the minimum
group size. The user can specify the minimum number of samples in which
the gene should be expressed. Essentially, a gene is kept only if it has
at least 1 count per million (CPM) in at least `min_group_size` samples.

## Usage

``` r
filter_low_expression_genes(raw_counts, min_group_size = 4)
```

## Arguments

- raw_counts:

  A data frame containing the raw counts of genes. The data frame should
  have samples as columns and genes as rows. All values should be
  non-negative integers.

- min_group_size:

  The minimum number of samples in which the gene should be expressed.
  Default is 4.

## Value

A data frame containing the filtered raw counts
