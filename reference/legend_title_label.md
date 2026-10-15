# Legend title for a colour scale

Legend title for a colour scale

## Usage

``` r
legend_title_label(
  legend.title,
  legend.title.face,
  scale_features,
  features,
  is_gene,
  suffix
)
```

## Arguments

- legend.title, legend.title.face:

  As in
  [`FeaturePlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/FeaturePlot_bioToolkit.md)

- scale_features:

  Features sharing the scale

- features:

  All features plotted

- is_gene:

  Logical vector named by feature

- suffix:

  Text appended to automatic titles, or NULL

## Value

NULL, a string, or a plotmath expression with gene names in italics
