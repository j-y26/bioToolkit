# Fill the holes of a logical mask

Fill the holes of a logical mask

## Usage

``` r
fill_holes(mask)
```

## Arguments

- mask:

  A logical matrix

## Value

The mask with every FALSE cell that cannot be reached from the border
without crossing a TRUE cell set to TRUE
