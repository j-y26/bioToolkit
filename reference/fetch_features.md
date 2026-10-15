# Fetch features for a fixed set of cells

Genes of the assay are fetched with its key prefix (e.g. "sct_Cdh5"), so
that they come from that assay only. Metadata columns and anything else
[`Seurat::FetchData()`](https://satijalab.github.io/seurat-object/reference/FetchData.html)
accepts (reduction dimensions, features already prefixed with a key) are
fetched as they are; a metadata column takes precedence over a gene of
the same name, as in FetchData().

## Usage

``` r
fetch_features(object, features, cells, slot, assay = NULL)
```

## Arguments

- object:

  A Seurat object

- features:

  Features to fetch

- cells:

  Cell names, in the order they are wanted

- slot:

  Layer to fetch expression from

- assay:

  The result of
  [`resolve_assay()`](https://j-y26.github.io/bioToolkit/reference/resolve_assay.md),
  or NULL to let FetchData() search the default assay

## Value

A data frame with one column per feature and one row per cell
