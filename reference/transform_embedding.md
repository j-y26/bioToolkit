# Tranforming a dimensionality reduction embedding, flipping and/or rotating

This function transforms a dimensionality reduction embedding by
flipping and/or rotating the embedding. This is useful for visualizing
the embedding in a consistent way, for example, when comparing two
embeddings or when the embedding is not in the desired orientation. This
function will only work for 2D embeddings. It will modify/overwrite the
original embedding matrix. However, the biological meaning of the
embedding is not changed. Only for visualization purposes.

## Usage

``` r
transform_embedding(
  obj,
  reduction = "umap",
  flip_horizontal = FALSE,
  flip_vertical = FALSE,
  rotate = 0
)
```

## Arguments

- obj:

  A Seurat object

- reduction:

  The name of the dimensionality reduction to be transformed (e.g.
  "pca", "umap", "tsne")

- flip_horizontal:

  A boolean indicating whether to flip the embedding horizontally
  (default = FALSE)

- flip_vertical:

  A boolean indicating whether to flip the embedding vertically (default
  = FALSE)

- rotate:

  A numeric value indicating the angle (in degrees) to rotate the
  embedding (default = 0)
