**bioToolkit** collects utility functions for transcriptomic and
epigenomic analysis in R: processing **single-cell** and **bulk** data,
summarising differential expression and gene-set enrichment, and turning
the results into **publication-quality figures**, built on ggplot2, that
look alike and combine cleanly.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("j-y26/bioToolkit")
```

Packages needed by only some functions (Seurat, edgeR, ggrepel, …) are
suggested, not required; a function that needs one says which.

## Guides

Each guide answers one question with public data, and shows every mode
of its functions with real output.

[**Figure essentials and
colours**](https://j-y26.github.io/bioToolkit/articles/figures.md)

Theme, fixed panel sizes, saving to PDF / PNG, legends as files, source
data, palettes ([all
palettes](https://j-y26.github.io/bioToolkit/articles/palettes.md)).

[`save_plot()`](https://j-y26.github.io/bioToolkit/reference/save_plot.md)
[`theme_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/theme_bioToolkit.md)
[`extract_legend()`](https://j-y26.github.io/bioToolkit/reference/extract_legend.md)
[`bt_palette()`](https://j-y26.github.io/bioToolkit/reference/bt_palette.md)

[**Single-cell:
embeddings**](https://j-y26.github.io/bioToolkit/articles/sc-embeddings.md)

Where cell types, conditions, genes and scores sit on a UMAP.

[`DimPlot_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/DimPlot_bioToolkit.md)
[`FeaturePlot_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/FeaturePlot_bioToolkit.md)
[`format_sc_plots()`](https://j-y26.github.io/bioToolkit/reference/format_sc_plots.md)

[**Single-cell: expression across
groups**](https://j-y26.github.io/bioToolkit/articles/sc-expression.md)

Marker genes and expression per cell type or condition.

[`top_markers()`](https://j-y26.github.io/bioToolkit/reference/top_markers.md)
[`DotPlot_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/DotPlot_bioToolkit.md)
[`VlnPlot_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/VlnPlot_bioToolkit.md)
[`DoHeatmap_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/DoHeatmap_bioToolkit.md)

[**Single-cell: comparing
conditions**](https://j-y26.github.io/bioToolkit/articles/sc-conditions.md)

Cell type composition, and where on the embedding cells are gained or
lost.

[`plot_composition_scatter()`](https://j-y26.github.io/bioToolkit/reference/plot_composition_scatter.md)
[`DensityDiffPlot_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/DensityDiffPlot_bioToolkit.md)

[**Bulk RNA-seq: samples and
structure**](https://j-y26.github.io/bioToolkit/articles/bulk-samples.md)

Filtering, PCA and what drives the variation between samples.

[`plot_pca()`](https://j-y26.github.io/bioToolkit/reference/plot_pca.md)
[`plot_variance_explained()`](https://j-y26.github.io/bioToolkit/reference/plot_variance_explained.md)
[`calc_expr_corr()`](https://j-y26.github.io/bioToolkit/reference/calc_expr_corr.md)

[**Differential expression: one
comparison**](https://j-y26.github.io/bioToolkit/articles/de-single.md)

edgeR, DESeq2, limma or Seurat results as volcano and MA plots, gene
lists and rankings.

[`standardize_de()`](https://j-y26.github.io/bioToolkit/reference/standardize_de.md)
[`plot_volcano()`](https://j-y26.github.io/bioToolkit/reference/plot_volcano.md)
[`plot_ma()`](https://j-y26.github.io/bioToolkit/reference/plot_ma.md)
[`filter_de_genes()`](https://j-y26.github.io/bioToolkit/reference/filter_de_genes.md)

[**Differential expression: many
comparisons**](https://j-y26.github.io/bioToolkit/articles/de-multi.md)

DE across cell types, treatments or assays, and whether two responses
agree.

[`plot_multi_volcano()`](https://j-y26.github.io/bioToolkit/reference/plot_multi_volcano.md)
[`plot_de_dots()`](https://j-y26.github.io/bioToolkit/reference/plot_de_dots.md)
[`plot_fc_quadrant()`](https://j-y26.github.io/bioToolkit/reference/plot_fc_quadrant.md)

[**Gene-set
enrichment**](https://j-y26.github.io/bioToolkit/articles/enrichment.md)

Over-representation and GSEA results from clusterProfiler, fgsea or
GSEA.

[`plot_enrichment_dots()`](https://j-y26.github.io/bioToolkit/reference/plot_enrichment_dots.md)
[`plot_enrichment_grid()`](https://j-y26.github.io/bioToolkit/reference/plot_enrichment_grid.md)
[`plot_gsea_rank()`](https://j-y26.github.io/bioToolkit/reference/plot_gsea_rank.md)

All functions, grouped the same way, are listed in the [function
index](https://j-y26.github.io/bioToolkit/reference/index.md).

## Conventions in 30 seconds

|  |  |
|----|----|
| `group.by` / `split.by` | `group.by`: the groups compared, marked by their colour; `split.by`: a higher-level grouping drawn as panels or labelled blocks, each marked by a bar. `group.by = "condition", split.by = "cell_type"` compares conditions within each cell type |
| `cols` | groups: a palette name or colours named by group; continuous scales: a palette name or colours from low to high. Diverging scales keep their centre at 0, each side reaching its full colour at its own extreme |
| sizes | fixed panels, the same in every figure: `panel.width` / `panel.height` for continuous axes, `cell.width` / `cell.height` (one row or column) for categories, in inches |
| look | [`theme_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/theme_bioToolkit.md): 11 / 9 / 12 pt text, 1 pt axes, bold legend titles; gene names in italics |
| background | transparent (`bg = "transparent"`), so figures sit on any slide or panel |
| saving | `save_plot(p, "fig.pdf")` sizes the file to fit: PDF (vector) or PNG (1000 dpi by default) |
| legends | `extract_legend(p)` takes a legend out, to save on its own |
| data | `export_source_data(p, "fig.xlsx")` writes the table behind a figure |
| DE tables | edgeR, DESeq2, limma and Seurat columns are recognised; others are named with `logfc.col`, `padj.col`, … |
