# Extract PCA loadings from PCA object

This function extracts the PCA loadings from a PCA object generated
using the `prcomp` function.

## Usage

``` r
extract_pc_loadings(pca_obj, pcs = 1:3, top_n = 20)
```

## Arguments

- pca_obj:

  A PCA object generated using the `prcomp` function.

- pcs:

  A vector of integers specifying the principal components to extract
  the loadings. Default is to extract PC 1-3.

- top_n:

  The number of top genes to extract for each principal component.
  Default is 20.

## Value

A list, where each element is a PC, and each element is a named vector
containing the top genes and their loadings.
