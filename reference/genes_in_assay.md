# Split genes into those of an assay and those missing from it

Missing genes are reported in a warning, naming those that match a gene
of the assay only in case (gene symbols differ in case between species:
human CD3E, mouse and rat Cd3e) and the assays that do hold the others.

## Usage

``` r
genes_in_assay(object, genes, assay)
```

## Arguments

- object:

  A Seurat object

- genes:

  Gene names

- assay:

  The result of
  [`resolve_assay()`](https://j-y26.github.io/bioToolkit/reference/resolve_assay.md)

## Value

The genes found in the assay, in the order given
