# bioToolkit 0.7.0

Phase C of the visualization overhaul: the dot-grid family and DE summaries.

## New

* `DotPlot_bioToolkit()`: marker dot plot (% of cells expressing = dot area,
  mean expression = colour), covering `Seurat::DotPlot()` and
  `SCpubr::do_DotPlot()`: `scale` (z-scored means, default TRUE),
  `average = "linear"` (Seurat: `log1p(mean(expm1(x)))`) or `"log"` (SCpubr:
  `mean(x)`), gene groups from a named list (labelled brackets), `split.by`
  (rows of a group kept together), clustering of genes and groups, `flip`,
  `pct.min` (SCpubr's `dot.min`) and `pct.max`. New: group markers in the
  groups' `DimPlot_bioToolkit()` colours (`group.marker`: a coloured dot by
  default, a colour strip, or none), unscaled means on "viridis_r" (yellow =
  low), gene names in italics, dot
  areas on a fixed 0-100% scale (comparable between figures), fixed cell
  sizes, and source data with the assay and layer used.
* Dot-grid plots (`DotPlot_bioToolkit()`, `plot_de_dots()`,
  `plot_enrichment_grid()`) have 1 pt x and y axes on the dots, with the
  group markers outside the axes and no grid lines through them;
  `theme_bioToolkit("grid")` has the same 1 pt axes.
* `top_markers()`: top marker genes per group from `FindAllMarkers()` output,
  as the named list `DotPlot_bioToolkit()` takes.
* `plot_de_dots()`: genes x comparisons, colour = log2 fold change, size =
  -log10 adjusted p-value; non-significant results faded, significant ones
  ringed, in a "Significance" legend that spells out both cutoffs
  (`significance = "fade"`, or "grey" / "hide");
  reads edgeR, DESeq2, limma and Seurat results through `standardize_de()`;
  gene groups, comparison groups and clustering. Replaces `de_dot_plot()`.
* `plot_multi_volcano()`: all DE genes of many comparisons, one column each,
  with a band at +-`logfc.cutoff`, counts of up / down genes, labels of the
  top genes per column; points spread by density (`layout = "sina"`), by
  p-value or at random; coloured by significance, direction or the number of
  comparisons sharing the gene; optional groups of comparisons.
* `standardize_enrichment()`: reads fgsea, clusterProfiler (GSEA and
  over-representation, objects or tables) and GSEA desktop results into one
  table; MSigDB names made readable ("TNFA signaling via NFKB"), GO
  descriptions capitalised ("Extracellular matrix organization").
* `plot_enrichment_dots()`: the top gene sets of one result, or of several
  (one panel each, e.g. up- and down-regulated genes, from a named list):
  gene ratio (or NES) on x, size = count, colour = -log10 padj. Each panel
  has its own single-hue palette and a tinted strip; panels stacked (shared
  x-axis), side by side or in a grid; shared colour / size limits by default;
  legends always in the same order (colour bar, then size); `term.replace`
  to shorten names. GSEA counts and ratios come from the leading edge.
* `plot_enrichment_grid()`: gene sets x comparisons (NES or -log10 padj),
  non-significant results faded.
* `plot_enrichment_bars()`: diverging NES bar chart of one comparison, names
  opposite their bars; faded (non-significant) bars explained in a legend.
* `plot_gsea_rank()`: running enrichment score, gene positions over a band
  coloured by the ranking metric, and the ranked metric (grey; `show.metric
  = FALSE` drops it) of one gene set, with NES and FDR; compact, touching
  panels. Takes a clusterProfiler `gseaResult`, or fgsea results with the
  ranked list and gene sets. The score is computed in the package and checked
  against the reported one.
* Palettes: single-hue sequential "seq-purple", "seq-blue", "seq-teal",
  "seq-green", "seq-orange", "seq-red", "seq-magenta" (from
  `grDevices::hcl.colors()`, trimmed to the same lightness range). Any
  palette name with "_r" appended is that palette reversed ("viridis_r",
  "npg_r", "seq-blue_r", and RColorBrewer names such as "BuPu_r");
  reversed palettes are not stored.

## Changed

* The sequential palette "bioToolkit" is removed: it was viridisLite's
  plasma reversed (and trimmed at its pale end). `FeaturePlot_bioToolkit()`
  now defaults to "plasma_r", plasma exactly reversed, so its lowest values
  are plasma's pale yellow (#F0F921) rather than the golden yellow the
  trimmed version ended on.

* `FeaturePlot_bioToolkit()` gains `assay`: genes are read from that assay
  only (the default assay when NULL, as before), the assay and layer are
  reported in a message and written to the source data, and a gene that is
  only in another assay is an error instead of being read from it silently.
* `save_plot()`: a legend taller (or wider) than its panel no longer adds
  blank space around single plots. The file is enlarged only as far as the
  legend actually spills beyond the plot (as already done for patchworks).
  Files of plots whose legend fits are unchanged.
* `fetch_features()` and `feature_is_gene()` moved to `R/seurat_data.R`, with
  the new assay helpers.
* `DOSE` added to Suggests (reading clusterProfiler result objects).

## Removed

* `de_dot_plot()`: use `plot_de_dots()`. Its arguments map as `de_results` ->
  `de` (with `comparison.col`), `feature_var` -> `feature.col`, `logFC_var` ->
  `logfc.col`, `p_val_var` -> `padj.col`, `logFC_cutoff` / `p_val_cutoff` ->
  `logfc.cutoff` / `padj.cutoff`, `invert` -> `flip`, `sig_only = TRUE` ->
  `significance = "hide"`, `group_by_features` -> `features` as a named list,
  `group_by_comparisons` -> `groups`, colours -> `cols`.


# bioToolkit 0.6.0

Phase B of the visualization overhaul: fixes and removals.

## New

* `extract_legend()`: takes the legend out of a plot, to save it as its own
  file with `save_plot()` (sized to fit, with the fonts of the output format)
  or place it in a layout. By default the legend is exactly as in the plot;
  `direction`, `title.position` and any `legend.*` theme setting change it.
  Hidden legends (`legend.position = "none"`) can be extracted too.
* `save_plot()` accepts legends from `extract_legend()`.

## Removed

* `annotation_sankey_plot()`. Annotation concordance is better shown as a
  contingency heatmap (e.g. with ComplexHeatmap). `plotly` is no longer
  suggested.
* `plot_top_pca_loadings()`. `extract_pc_loadings()` stays; a PCA score plot
  (`plot_pca()`) is planned.
* `plot_volcano2()`: merged into `plot_volcano()` (`annotate.counts = "arrow"`).

## Changed

* The sequential palette "bioToolkit" is removed: it was viridisLite's
  plasma reversed (and trimmed at its pale end). `FeaturePlot_bioToolkit()`
  now defaults to "plasma_r", plasma exactly reversed, so its lowest values
  are plasma's pale yellow (#F0F921) rather than the golden yellow the
  trimmed version ended on.

* `plot_volcano()` is rewritten (breaking change: new arguments):
  * reads edgeR, DESeq2, limma and Seurat results through `standardize_de()`
    (`method`, `feature.col`, `logfc.col`, `p.col`, `padj.col` replace
    `de_method`, `x`, `y`, `label`);
  * the y-axis is the -log10 adjusted p-value, the same value that decides
    significance (`padj.cutoff`, `logfc.cutoff` replace `p_cutoff`,
    `fc_cutoff`);
  * labels are chosen explicitly with `label = "top"` (the `label.n` most
    significant per direction, plus any `label.features`), `"selected"` (only
    `label.features`) or `"none"`. The old behaviour of silently replacing
    requested genes with the top 10 per direction is gone, and labels are no
    longer dropped when crowded. Labels are bold, and italic for genes
    (`italic.genes`);
  * `annotate.counts = "corner"` (old `plot_volcano()`), `"arrow"` (old
    `plot_volcano2()`) or `"none"`;
  * fixed panel size (`panel.width`, `panel.height`; square by default:
    `panel.height = panel.width`, so `panel.width` alone resizes the plot),
    `theme_bioToolkit()`,
    transparent background; points are rasterized above 5,000 features
    (`raster`, `raster.dpi = 1000`), labels and text stay vector;
  * non-significant points are drawn first, so they never hide DE points;
  * p-values of 0 no longer stretch the axis to 1e-300;
  * `cols = c(up = , down = , ns = )` replaces `up_color`, `down_color`,
    `non_de_color`; `up.label`, `down.label`, `ns.label` replace `up_label`,
    `down_label`, `non_de_label`. Text sizes follow the theme.
* `de_dot_plot()` fixes: non-significant dots are now grey (`nsig_color`,
  default `"grey70"` as documented); logFC values beyond `col_limits_logFC`
  take the colour of the nearest limit instead of the non-significant grey;
  the default colour limits are symmetric and computed from the plotted
  significant dots (`col_limits_logFC = NULL`); p-values of 0 no longer give
  infinite dot sizes; features are matched case-sensitively (gene-symbol case
  differs between species); the size legend names the p-value column; the
  group strips of the inverted plot follow `x_group_label_position`; the
  "removed due to non-significance" warning only appears with `sig_only = TRUE`.
* `calc_expr_corr()`: its volcano works again (it passed a vector as `label`),
  is saved at the volcano's own fixed size, and `base_dir` is created even when
  `plot_model = FALSE`. Default `dpi` is now 1000.
* Theme size presets are removed: the "nature", "cell" and "slides" presets,
  the option `bioToolkit.preset` and the `preset` argument of
  `theme_bioToolkit()`. One set of sizes applies everywhere: 11 pt text (9 pt
  tick labels and legend text, 12 pt titles and strips), 1 pt axis lines and
  ticks (were 0.5 pt), 0.5 pt other lines. `theme_bioToolkit(base.size = )`
  still scales the text.
* `DimPlot_bioToolkit()` (`label.size`, `legend.pt.size`) and
  `FeaturePlot_bioToolkit()` (`label.size`, `legend.length`, `legend.width`)
  show their defaults again (3.5, 3.5; 3.5, 1.5, 0.15) instead of `NULL`;
  the plots are unchanged.
* Function signatures list one argument per line (no change in behaviour).

# bioToolkit 0.5.0

Phase A of the visualization overhaul: the foundations every plot builds on.

## New

* `theme_bioToolkit()`: one theme for all plots, in three layouts ("axis",
  "grid", "embedding"), with size presets for the printed figure ("default",
  "nature", "cell", "slides"). `options(bioToolkit.preset = "nature")` switches
  every plot of the package at once.
* `italic_labels()`: gene names in italics as scale labels.
* `panel_size_theme()`: fixed panel sizes in inches for any ggplot.
* `save_plot()`: saves a plot to a file sized to fit it (replaces
  `save_dimplot()`); also saves plots without fixed panels at a given `width`
  and `height`.
* Palettes: `bt_palette()` returns the colours of a palette as a plain hex
  vector, to pick, reorder and name by index
  (`bt_palette("scanpy-20")[c(1, 3, 2)]`); `show_palettes()` shows each colour
  with its index; `list_palettes()` lists names, types and sources;
  `check_palette()` flags colours that are hard to tell apart with colour
  vision deficiency. Ships the scanpy palettes, colour-blind safe palettes
  (Okabe-Ito, Paul Tol), journal styles (NPG, AAAS, NEJM, Lancet, JAMA, JCO),
  viridis and Tol sequential palettes, and diverging palettes.
* `standardize_de()`: reads edgeR, DESeq2, limma and Seurat results (or a list
  of them, one per comparison) into standard columns.
* `export_source_data()`: writes the data behind a plot, or each plot of a
  patchwork, to CSV or XLSX.

## Changed

* The sequential palette "bioToolkit" is removed: it was viridisLite's
  plasma reversed (and trimmed at its pale end). `FeaturePlot_bioToolkit()`
  now defaults to "plasma_r", plasma exactly reversed, so its lowest values
  are plasma's pale yellow (#F0F921) rather than the golden yellow the
  trimmed version ended on.

* Plots have a transparent background by default (`bg = "transparent"`), so
  that they can be placed on any background; use `bg = "white"` for an opaque
  one. `save_plot()` keeps it in PDF, PNG, TIFF and SVG, including for
  patchworks, and uses white for JPEG and BMP.
* The default colours of groups are now "scanpy-20" for up to 20 groups and
  "scanpy-102" beyond (before: the ggplot2 hue palette). Every level of a
  grouping variable gets a colour, including levels with no plotted cell, so
  that a plot of a subset is coloured like a plot of the whole object.
* `cols` of `DimPlot_bioToolkit()` accepts palette names; `cols` of
  `FeaturePlot_bioToolkit()` accepts palette names or colours from low to
  high (`c(low, high)`, `c(low, mid, high)`). Colours given as a function and
  the viridis letters "A" to "H" are no longer accepted.
* `DimPlot_bioToolkit()` and `FeaturePlot_bioToolkit()` take their text sizes,
  legend keys, colour bar and arrows from the theme preset: `legend.text.size`,
  `strip.text.size`, `label.size`, `legend.pt.size`, `legend.length` and
  `legend.width` now default to `NULL` (from the preset). At the default preset
  the plots are unchanged.
* `grouped_bar_plot()` and `peak_plot_frip()` take their colours from the
  package palettes instead of Seurat (default `"polychrome-36"`, the same colours as
  Seurat's "polychrome"); Seurat's "glasbey", "stepped" and "parade" are no
  longer available.
* `calc_expr_corr()` saves its volcano plot with `save_plot()`.
* Requires R >= 4.1.0.

## Removed

* `save_dimplot()` (renamed `save_plot()`) and its deprecated alias
  `save_umap()`.

## Dependencies

* Suggests: `writexl` (XLSX export); `png` and `jpeg` (tests only).
