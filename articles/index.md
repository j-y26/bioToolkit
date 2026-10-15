# Articles

### Figures

- [Figure essentials and
  colours](https://j-y26.github.io/bioToolkit/articles/figures.md):

  The conventions every bioToolkit plot shares: theme, fixed panel
  sizes, saving to PDF and PNG, legends, source data, gene names and
  colours.

- [Colour
  palettes](https://j-y26.github.io/bioToolkit/articles/palettes.md):

  Every palette of bioToolkit with the index and hex code of each
  colour, its source, and the plots that use it by default.

### Single-cell

- [Single-cell:
  embeddings](https://j-y26.github.io/bioToolkit/articles/sc-embeddings.md):

  Where cell types, conditions, genes and scores sit on a UMAP:
  DimPlot_bioToolkit(), FeaturePlot_bioToolkit(), orienting an embedding
  and restyling plots made by Seurat.

- [Single-cell: expression across
  groups](https://j-y26.github.io/bioToolkit/articles/sc-expression.md):

  Marker genes and expression per cell type or condition: top_markers(),
  DotPlot_bioToolkit(), VlnPlot_bioToolkit(), RidgePlot_bioToolkit() and
  DoHeatmap_bioToolkit().

- [Single-cell: comparing
  conditions](https://j-y26.github.io/bioToolkit/articles/sc-conditions.md):

  Does a condition change which cells are present, or where they sit on
  the embedding? summarize_group_meta(), get_meta_pivot(),
  plot_composition_scatter() and DensityDiffPlot_bioToolkit().

### Bulk and differential expression

- [Bulk RNA-seq: samples and
  structure](https://j-y26.github.io/bioToolkit/articles/bulk-samples.md):

  From a count matrix to the structure of the samples: sample
  annotation, low-expression filtering, PCA, what drives the variation,
  group summaries and genes co-expressed with a gene of interest.

- [Differential expression: one
  comparison](https://j-y26.github.io/bioToolkit/articles/de-single.md):

  Read a DE table from edgeR, DESeq2, limma or Seurat; volcano and MA
  plots; up- and down-regulated gene lists and ranked lists for GSEA.

- [Differential expression: many
  comparisons](https://j-y26.github.io/bioToolkit/articles/de-multi.md):

  Compare DE results across cell types, treatments, time points or
  assays: plot_multi_volcano(), plot_de_dots() and plot_fc_quadrant().

### Gene sets

- [Gene-set
  enrichment](https://j-y26.github.io/bioToolkit/articles/enrichment.md):

  Plot over-representation and GSEA results from clusterProfiler, fgsea
  or the GSEA application: standardize_enrichment(),
  plot_enrichment_dots(), plot_enrichment_bars(),
  plot_enrichment_grid(), plot_gsea_rank() and format_gene_sets().
