# Resolve the assay and layer expression is read from

Resolve the assay and layer expression is read from

## Usage

``` r
resolve_assay(object, assay, slot, quiet = FALSE)
```

## Arguments

- object:

  A Seurat object

- assay:

  Assay name, or NULL for the default assay

- slot:

  Layer name (e.g. "data", "counts")

- quiet:

  Whether to skip the message naming the assay and layer

## Value

A list with the \`assay\` name, its \`key\` (e.g. "sct\_"), the
\`layer\`, its \`features\` (row names) and whether the assay was the
\`default\`
