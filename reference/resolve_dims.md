# Resolve the dimensions to plot and their arrow labels

Before \`dims\` selected the dimensions, it was passed through \`...\`
to the arrows as their labels. A character \`dims\` is therefore still
taken as the arrow labels of the first two dimensions.

## Usage

``` r
resolve_dims(dims, reduction)
```

## Arguments

- dims:

  Two dimension indices, or two arrow labels

- reduction:

  Name of the reduction

## Value

A list with the integer \`dims\` and their arrow \`labels\`
