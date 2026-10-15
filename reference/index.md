# Package index

## Figure essentials

Shared by every plot: the theme, fixed panel sizes, saving to PDF and
PNG, legends as separate files, the data behind a plot, italic gene
names. Guide: [Figure essentials and
colours](https://j-y26.github.io/bioToolkit/articles/figures.md).

- [`theme_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/theme_bioToolkit.md)
  : The bioToolkit ggplot2 theme
- [`panel_size_theme()`](https://j-y26.github.io/bioToolkit/reference/panel_size_theme.md)
  : Theme giving panels a fixed absolute size
- [`save_plot()`](https://j-y26.github.io/bioToolkit/reference/save_plot.md)
  : Save a plot to a file sized to fit it
- [`extract_legend()`](https://j-y26.github.io/bioToolkit/reference/extract_legend.md)
  [`print(`*`<bt_legend>`*`)`](https://j-y26.github.io/bioToolkit/reference/extract_legend.md)
  : Extract the legend of a plot
- [`export_source_data()`](https://j-y26.github.io/bioToolkit/reference/export_source_data.md)
  : Export the data behind a plot
- [`italic_labels()`](https://j-y26.github.io/bioToolkit/reference/italic_labels.md)
  : Labels with gene names in italics

## Colours

Browse, pick and check palettes; every plot takes one `cols` argument.
All palettes: [Colour
palettes](https://j-y26.github.io/bioToolkit/articles/palettes.md).
Guide: [Figure essentials and
colours](https://j-y26.github.io/bioToolkit/articles/figures.html#colours).

- [`bt_palette()`](https://j-y26.github.io/bioToolkit/reference/bt_palette.md)
  : Colours of a palette
- [`list_palettes()`](https://j-y26.github.io/bioToolkit/reference/list_palettes.md)
  : List the palettes of the package
- [`show_palettes()`](https://j-y26.github.io/bioToolkit/reference/show_palettes.md)
  : Show palettes with the index of each colour
- [`check_palette()`](https://j-y26.github.io/bioToolkit/reference/check_palette.md)
  : Check that colours stay distinct with colour vision deficiency
- [`mix_colors()`](https://j-y26.github.io/bioToolkit/reference/mix_colors.md)
  : Mix colors based on weights of each color
- [`sort_colors_by_hue()`](https://j-y26.github.io/bioToolkit/reference/sort_colors_by_hue.md)
  : Sort the Colors by Hue
- [`ggplot_colors()`](https://j-y26.github.io/bioToolkit/reference/ggplot_colors.md)
  : Getting default ggplot2 colors

## Single-cell: embeddings

Cells, groups and genes on a UMAP / t-SNE; orienting embeddings;
restyling plots made by Seurat. Guide:
[Embeddings](https://j-y26.github.io/bioToolkit/articles/sc-embeddings.md).

- [`DimPlot_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/DimPlot_bioToolkit.md)
  : Publication-quality dimensionality reduction plot
- [`FeaturePlot_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/FeaturePlot_bioToolkit.md)
  : Publication-quality feature plot
- [`transform_embedding()`](https://j-y26.github.io/bioToolkit/reference/transform_embedding.md)
  : Tranforming a dimensionality reduction embedding, flipping and/or
  rotating
- [`format_sc_plots()`](https://j-y26.github.io/bioToolkit/reference/format_sc_plots.md)
  : Format a set of embedding plots (from Seurat, scCustomize, SCpubr
  etc.)
- [`manipulate_sc_plot()`](https://j-y26.github.io/bioToolkit/reference/manipulate_sc_plot.md)
  : Restyle a single ggplot made by Seurat
- [`add_dimension_arrows()`](https://j-y26.github.io/bioToolkit/reference/add_dimension_arrows.md)
  : Add dimension arrows

## Single-cell: expression across groups

Marker genes and expression per cluster, cell type or condition. Guide:
[Expression across
groups](https://j-y26.github.io/bioToolkit/articles/sc-expression.md).

- [`top_markers()`](https://j-y26.github.io/bioToolkit/reference/top_markers.md)
  : Top marker genes per group
- [`DotPlot_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/DotPlot_bioToolkit.md)
  : Marker dot plot
- [`VlnPlot_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/VlnPlot_bioToolkit.md)
  : Violin plot of genes or per-cell values
- [`RidgePlot_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/RidgePlot_bioToolkit.md)
  : Ridge plot of genes or per-cell values
- [`DoHeatmap_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/DoHeatmap_bioToolkit.md)
  : Heatmap of genes across cells or groups

## Single-cell: comparing conditions

Cell type composition and where cells gain or lose density between
conditions. Guide: [Comparing
conditions](https://j-y26.github.io/bioToolkit/articles/sc-conditions.md).

- [`summarize_group_meta()`](https://j-y26.github.io/bioToolkit/reference/summarize_group_meta.md)
  : Count and proportion of cells per group and condition
- [`plot_composition_scatter()`](https://j-y26.github.io/bioToolkit/reference/plot_composition_scatter.md)
  : Composition scatter: group proportions in two conditions
- [`DensityDiffPlot_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/DensityDiffPlot_bioToolkit.md)
  : Difference in cell density between conditions on an embedding
- [`get_meta_pivot()`](https://j-y26.github.io/bioToolkit/reference/get_meta_pivot.md)
  : Generate a pivot table from cell metadata

## Bulk RNA-seq: samples and structure

Sample annotation, low-expression filtering, PCA, what drives variation,
group summaries and co-expression with a gene. Guide: [Bulk
RNA-seq](https://j-y26.github.io/bioToolkit/articles/bulk-samples.md).

- [`keywords_to_info()`](https://j-y26.github.io/bioToolkit/reference/keywords_to_info.md)
  : Extracting sample info based on regex keywords
- [`filter_low_expression_genes()`](https://j-y26.github.io/bioToolkit/reference/filter_low_expression_genes.md)
  : Filter out genes with low expression based on minimum group size
- [`plot_pca()`](https://j-y26.github.io/bioToolkit/reference/plot_pca.md)
  : PCA score plot of samples
- [`extract_pc_loadings()`](https://j-y26.github.io/bioToolkit/reference/extract_pc_loadings.md)
  : Extract PCA loadings from PCA object
- [`plot_variance_explained()`](https://j-y26.github.io/bioToolkit/reference/plot_variance_explained.md)
  : Variance explained by sample covariates
- [`summarize_group_expression()`](https://j-y26.github.io/bioToolkit/reference/summarize_group_expression.md)
  : Summarize group expression values for each gene
- [`calc_expr_corr()`](https://j-y26.github.io/bioToolkit/reference/calc_expr_corr.md)
  : Calculate the correlation between gene expression and another
  variable using a linear model

## Differential expression: one comparison

Read a DE table from edgeR, DESeq2, limma or Seurat; volcano and MA
plots; up / down gene lists and ranked lists for GSEA. Guide: [DE: one
comparison](https://j-y26.github.io/bioToolkit/articles/de-single.md).

- [`standardize_de()`](https://j-y26.github.io/bioToolkit/reference/standardize_de.md)
  : Standardize differential expression results
- [`plot_volcano()`](https://j-y26.github.io/bioToolkit/reference/plot_volcano.md)
  : Volcano plot of differential expression results
- [`plot_ma()`](https://j-y26.github.io/bioToolkit/reference/plot_ma.md)
  : MA plot of differential expression results
- [`filter_de_genes()`](https://j-y26.github.io/bioToolkit/reference/filter_de_genes.md)
  : Split DE results into up- and down-regulated genes
- [`rank_de_genes()`](https://j-y26.github.io/bioToolkit/reference/rank_de_genes.md)
  : Calculate ranks of genes for differential expression analysis
- [`generate_de_rank_file()`](https://j-y26.github.io/bioToolkit/reference/generate_de_rank_file.md)
  : Generate rank file for GSEA

## Differential expression: many comparisons

DE across cell types, time points, treatments or assays. Guide: [DE:
many
comparisons](https://j-y26.github.io/bioToolkit/articles/de-multi.md).

- [`plot_multi_volcano()`](https://j-y26.github.io/bioToolkit/reference/plot_multi_volcano.md)
  : Volcano-style summary of many DE comparisons
- [`plot_de_dots()`](https://j-y26.github.io/bioToolkit/reference/plot_de_dots.md)
  : Dot plot of differential expression across comparisons
- [`plot_fc_quadrant()`](https://j-y26.github.io/bioToolkit/reference/plot_fc_quadrant.md)
  : Quadrant plot of the log fold changes of two DE comparisons

## Gene-set enrichment

Over-representation and GSEA results from clusterProfiler, fgsea or
GSEA: dot plots, NES bars, grids across comparisons, running scores.
Guide: [Gene-set
enrichment](https://j-y26.github.io/bioToolkit/articles/enrichment.md).

- [`standardize_enrichment()`](https://j-y26.github.io/bioToolkit/reference/standardize_enrichment.md)
  : Standardize pathway enrichment results
- [`plot_enrichment_dots()`](https://j-y26.github.io/bioToolkit/reference/plot_enrichment_dots.md)
  : Dot plot of the top gene sets of enrichment results
- [`plot_enrichment_bars()`](https://j-y26.github.io/bioToolkit/reference/plot_enrichment_bars.md)
  : Diverging bar chart of normalized enrichment scores
- [`plot_enrichment_grid()`](https://j-y26.github.io/bioToolkit/reference/plot_enrichment_grid.md)
  : Grid of pathway enrichment across comparisons
- [`plot_gsea_rank()`](https://j-y26.github.io/bioToolkit/reference/plot_gsea_rank.md)
  : GSEA running enrichment score plot
- [`format_gene_sets()`](https://j-y26.github.io/bioToolkit/reference/format_gene_sets.md)
  : Convert gene sets between formats
