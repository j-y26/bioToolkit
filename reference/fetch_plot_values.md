# Fetch any numeric features of the cells: genes, metadata or embeddings

Used by plots that accept a mix of genes and per-cell values (violin,
ridge and heatmap plots). From a Seurat object, metadata columns (QC
metrics, module scores) are read as they are; genes are read from one
assay, keyed so that they never come from another one
([`fetch_features()`](https://j-y26.github.io/bioToolkit/reference/fetch_features.md)),
and the assay and layer are reported only when at least one feature is a
gene. Genes missing from the assay are left out with a warning that
names near matches in case
([`genes_in_assay()`](https://j-y26.github.io/bioToolkit/reference/genes_in_assay.md)).
From a data frame (e.g. a table of QC metrics), features are its
columns.

## Usage

``` r
fetch_plot_values(object, features, cells, assay, slot)
```

## Arguments

- object:

  A Seurat object, or a data frame with one row per cell

- features:

  Feature names

- cells:

  Cell names (Seurat object) or row names / indices (data frame), in the
  order wanted

- assay:

  Assay name, or NULL for the default assay

- slot:

  Layer to read genes from

## Value

A list with \`data\` (a data frame, one column per feature found and one
row per cell), \`is_gene\` (logical, named by feature) and \`assay\`
(the result of
[`resolve_assay()`](https://j-y26.github.io/bioToolkit/reference/resolve_assay.md),
or NULL when no feature is a gene)
