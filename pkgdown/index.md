<div class="bt-hero">

**bioToolkit** collects utility functions for transcriptomic and epigenomic
analysis in R: processing **single-cell** and **bulk** data, summarising
differential expression and gene-set enrichment, and turning the results into
**publication-quality figures**, built on ggplot2, that look alike and combine
cleanly.

</div>

## Installation

```r
# install.packages("remotes")
remotes::install_github("j-y26/bioToolkit")
```

Packages needed by only some functions (Seurat, edgeR, ggrepel, ...) are
suggested, not required; a function that needs one says which.

## Guides

Each guide answers one question with public data, and shows every mode of its
functions with real output.

::: {.bt-cards}

::: {.bt-card}
[**Figure essentials and colours**](articles/figures.html)

Theme, fixed panel sizes, saving to PDF / PNG, legends as files, source data,
palettes ([all palettes](articles/palettes.html)).

`save_plot()` `theme_bioToolkit()` `extract_legend()` `bt_palette()`
:::

::: {.bt-card}
[**Single-cell: embeddings**](articles/sc-embeddings.html)

Where cell types, conditions, genes and scores sit on a UMAP.

`DimPlot_bioToolkit()` `FeaturePlot_bioToolkit()` `format_sc_plots()`
:::

::: {.bt-card}
[**Single-cell: expression across groups**](articles/sc-expression.html)

Marker genes and expression per cell type or condition.

`top_markers()` `DotPlot_bioToolkit()` `VlnPlot_bioToolkit()` `DoHeatmap_bioToolkit()`
:::

::: {.bt-card}
[**Single-cell: comparing conditions**](articles/sc-conditions.html)

Cell type composition, and where on the embedding cells are gained or lost.

`plot_composition_scatter()` `DensityDiffPlot_bioToolkit()`
:::

::: {.bt-card}
[**Bulk RNA-seq: samples and structure**](articles/bulk-samples.html)

Filtering, PCA and what drives the variation between samples.

`plot_pca()` `plot_variance_explained()` `calc_expr_corr()`
:::

::: {.bt-card}
[**Differential expression: one comparison**](articles/de-single.html)

edgeR, DESeq2, limma or Seurat results as volcano and MA plots, gene lists
and rankings.

`standardize_de()` `plot_volcano()` `plot_ma()` `filter_de_genes()`
:::

::: {.bt-card}
[**Differential expression: many comparisons**](articles/de-multi.html)

DE across cell types, treatments or assays, and whether two responses agree.

`plot_multi_volcano()` `plot_de_dots()` `plot_fc_quadrant()`
:::

::: {.bt-card}
[**Gene-set enrichment**](articles/enrichment.html)

Over-representation and GSEA results from clusterProfiler, fgsea or GSEA.

`plot_enrichment_dots()` `plot_enrichment_grid()` `plot_gsea_rank()`
:::

:::

All functions, grouped the same way, are listed in the
[function index](reference/index.html).

## Conventions in 30 seconds

| | |
|---|---|
| `group.by` / `split.by` | `group.by`: the groups compared, marked by their colour; `split.by`: a higher-level grouping drawn as panels or labelled blocks, each marked by a bar. `group.by = "condition", split.by = "cell_type"` compares conditions within each cell type |
| `cols` | groups: a palette name or colours named by group; continuous scales: a palette name or colours from low to high. Diverging scales keep their centre at 0, each side reaching its full colour at its own extreme |
| sizes | fixed panels, the same in every figure: `panel.width` / `panel.height` for continuous axes, `cell.width` / `cell.height` (one row or column) for categories, in inches |
| look | `theme_bioToolkit()`: 11 / 9 / 12 pt text, 1 pt axes, bold legend titles; gene names in italics |
| background | transparent (`bg = "transparent"`), so figures sit on any slide or panel |
| saving | `save_plot(p, "fig.pdf")` sizes the file to fit: PDF (vector) or PNG (1000 dpi by default) |
| legends | `extract_legend(p)` takes a legend out, to save on its own |
| data | `export_source_data(p, "fig.xlsx")` writes the table behind a figure |
| DE tables | edgeR, DESeq2, limma and Seurat columns are recognised; others are named with `logfc.col`, `padj.col`, ... |
