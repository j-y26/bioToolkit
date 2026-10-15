# Convert gene sets between formats

Converts gene sets between a GMT file, a named list of gene vectors (the
input of fgsea and of
[`Seurat::AddModuleScore()`](https://satijalab.org/seurat/reference/AddModuleScore.html))
and a GSEABase \`GeneSetCollection\`.

## Usage

``` r
format_gene_sets(gsc, to = c("GeneSetCollection", "list", "gmt"), file = NULL)
```

## Arguments

- gsc:

  A GeneSetCollection object or a list of gene sets or a path to a GMT
  file.

- to:

  A character string indicating the format to convert to. Available
  options are "GeneSetCollection", "list", and "gmt".

- file:

  For \`to = "gmt"\`, a file to write the GMT lines to. Default is NULL
  (only returned)

## Value

A GeneSetCollection object, a named list of gene sets, or, for "gmt",
the GMT lines (one per set: name, an empty description, then the genes,
separated by tabs), invisibly when written to \`file\`

## Examples

``` r
sets <- list(SET_A = c("CD3E", "CD3D"), SET_B = c("MS4A1", "CD79A", "CD19"))
format_gene_sets(sets, to = "gmt")
#> [1] "SET_A\t\tCD3E\tCD3D"         "SET_B\t\tMS4A1\tCD79A\tCD19"
```
