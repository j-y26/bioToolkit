# A Seurat object, or a table as a data frame

Anything that is not a Seurat object is converted with
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) (data
frames, tibbles, data.tables, Bioconductor `DataFrame`s, `GRanges`).
Tables without row names of their own (automatic row names, or a
`GRanges` with duplicated names) are numbered, and flagged with the
attribute \`row_numbers\`.

## Usage

``` r
as_plot_table(object)
```

## Arguments

- object:

  A Seurat object or a table

## Value

The Seurat object, or a data frame with unique row names
