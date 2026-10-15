# Order the rows of a matrix by hierarchical clustering, within groups

Rows are clustered with Ward's method (\`ward.D2\`) on Euclidean
distances, separately within each group, and groups keep their order.
Missing values count as 0. Dendrograms are not drawn; only the order is
used.

## Usage

``` r
cluster_order(mat, groups = NULL)
```

## Arguments

- mat:

  A numeric matrix with row names

- groups:

  A factor along the rows, or NULL for one group

## Value

The row names, in clustered order
