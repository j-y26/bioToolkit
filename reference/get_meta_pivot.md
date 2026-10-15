# Generate a pivot table from cell metadata

This function generates a pivot table from cell metadata and returns a
data frame with the two cell metadata as rows and columns, with the
values being either the count of cells or the percentage of cells (by
rows or columns or total).

## Usage

``` r
get_meta_pivot(
  row_metadata,
  col_metadata,
  type = c("percent", "count"),
  margin = c("row", "column", "total")
)
```

## Arguments

- row_metadata:

  A vector of cell metadata to be used as rows, often a metadata that
  represents cell clusters or cell types

- col_metadata:

  A vector of cell metadata to be used as columns, must be the same
  length as row_metadata

- type:

  The type of pivot table to generate, either "percent" or "count"

- margin:

  The margin to calculate the percentage, either "row", "column" or
  "total", used only when type = "percent"

## Value

a data frame with the pivot table
