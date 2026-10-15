# Resolve the arrow arguments passed through \`...\` of an embedding plot

[`DimPlot_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/DimPlot_bioToolkit.md),
[`FeaturePlot_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/FeaturePlot_bioToolkit.md)
and
[`DensityDiffPlot_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/DensityDiffPlot_bioToolkit.md)
pass \`...\` to the dimension arrows. The arguments take the dot.case
names of
[`add_dimension_arrows`](https://j-y26.github.io/bioToolkit/reference/add_dimension_arrows.md)
(since 0.10.0); the old snake_case names are accepted with a deprecation
warning. \`arrow.scale\` (new in 1.0.0) has no old name. They are
returned under the snake_case names of the internal
[`dimension_arrow_grob`](https://j-y26.github.io/bioToolkit/reference/dimension_arrow_grob.md).

## Usage

``` r
resolve_arrow_args(dots, fn)
```

## Arguments

- dots:

  \`list(...)\` of the plot, without \`rasterize\`

- fn:

  Name of the plot function, for messages

## Value

A named list of arguments for
[`dimension_arrow_grob()`](https://j-y26.github.io/bioToolkit/reference/dimension_arrow_grob.md)
