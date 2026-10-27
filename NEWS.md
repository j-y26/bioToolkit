# bioToolkit 1.1.2

## Saving

* `save_plot()`: the `minus` argument (new in 1.1.0) is removed. Negative
  numbers on axes and colour bars keep ggplot2's hyphen ("-2"), as in PNGs,
  instead of the wider minus sign (U+2212).
* PDFs are still written with `cairo_pdf()`, so a hyphen stays a hyphen, but
  with cairo's default "sans" font, resolved by the system, rather than the
  font looked up with systemfonts (removed from Suggests).

# bioToolkit 1.1.1

Change default volcano up and down colors to brighter ones.

# bioToolkit 1.1.0

## Saving

* `save_plot()` writes PDFs with `cairo_pdf()` when R has cairo support.
  Base `pdf()` drew every hyphen with the wider minus glyph (its ISOLatin1
  encoding on Linux / macOS), could not draw Unicode text and did not embed
  its fonts. The PDF's default font is the one ragg uses for PNGs (found with
  systemfonts, new in Suggests), so PDF and PNG files of a plot use the same
  font and have nearly the same size. Without cairo, or with `device` given,
  `ggsave()`'s device is used as before.
* New `minus` argument (default TRUE): negative numbers in the labels of
  continuous and binned scales (axes, colour bars) are written with a true
  minus sign (U+2212). Text typed in titles, labels or gene names is not
  changed.

# bioToolkit 1.0.0

First stable release: plot revisions after the first round of the
documentation website, scaled dimension arrows, a palette page.

## Dimension arrows scale with the panel

* `DimPlot_bioToolkit()`, `FeaturePlot_bioToolkit()`,
  `DensityDiffPlot_bioToolkit()`, `format_sc_plots()`,
  `add_dimension_arrows()`: the dimension arrows keep a fixed proportion to
  the panel. Their sizes are those of a 3.5-inch panel and are scaled with
  the panel's shorter side when the plot is drawn: lengths, line width and
  label size. Small panels get small arrows that stay clear of the cells near
  the left and bottom edges (they used to run into them); large panels get
  larger arrows. New `arrow.scale` (default 1; in `...` of the embedding
  plots) makes them larger or smaller at any size. The geometry arguments
  (`arrow.length`, `text.font.size`, ...) are now sizes on a 3.5-inch panel.
  Plots without a fixed panel size (e.g. Seurat plots passed to
  `format_sc_plots()`) are scaled from the size they are drawn at.
* The room left below and left of embedding panels for the arrow labels
  scales with the arrows (0.6 cm at 3.5 inches).

## Group annotations: bars

Higher-level groups of rows or columns are drawn the same way in every plot:
a bar along the group's cells, 3.5 pt thick by default (`group.bar.width`),
set a little off the grid, with the group's name in bold beyond it. Bars are
black by default, or coloured by group (`*.cols`, the same values as `cols`).
By default they sit opposite the axis labels; `*.side = "labels"` moves them
beyond the labels. New shared code in `R/group_bars.R`.

* `DotPlot_bioToolkit()`: `feature.group.cols`, `feature.group.side`,
  `split.by.cols`, `split.by.side`, `group.bar.width`. Gene groups are now on
  top (right with `flip = TRUE`) with a bar instead of a thin bracket line;
  `split.by` blocks on the right (top when flipped) instead of the left.
* `VlnPlot_bioToolkit()`: `split.by.cols`, `split.by.side`,
  `group.bar.width`; flipped `split.by` blocks are now on the right.
* `DoHeatmap_bioToolkit()`: `feature.group.cols`, `split.by.cols`,
  `group.bar.width` (blocks of features opposite the feature names, blocks of
  columns on top).
* `plot_de_dots()`: comparisons are marked by colour like the groups of the
  dot plot (`comparison.marker` = "dot" (default), "strip" or "none";
  `comparison.cols`); `groups` and gene groups get bars (`groups.cols`,
  `groups.side`, `feature.group.cols`, `feature.group.side`,
  `group.bar.width`).
* The dot-grid size legend ("% expressing", "-log10 padj") is black.

## Changes

* `plot_multi_volcano()`: the coloured band at 0 is replaced by a strip of
  fixed thickness above each comparison's column with its name above it
  (`comparison.cols`, `comparison.labels`; formerly `band.cols`,
  `band.labels`, accepted with a deprecation warning), and dashed lines at the
  cutoffs. Non-significant genes are drawn in grey by default (`show.ns =
  TRUE`). `groups` get bars like `split.by` (`groups.cols`, `groups.side`).
  New `flip` lists the comparisons as rows (down left, up right; new
  `cell.height`, `panel.width`). Grouped comparisons are drawn together.
* Diverging colour scales (dot plot z-scores, heatmaps, DE dots, enrichment
  grid NES, multi-volcano, density difference, diverging feature plots) keep
  their centre colour at 0 but are no longer forced symmetric: each end is
  the cutoff, or the data's extreme on that side when closer, and each side
  is stretched to its own end (`diverging_limits()`, `rescale_diverging()`).
  Fixed `logfc.limits` given by the user are kept.
* Legend titles are bold in every plot (`theme_bioToolkit()`, also plotmath
  titles); `manipulate_sc_plot()` / `format_sc_plots()` default to
  `legend.title.bold = TRUE`.
* `DimPlot_bioToolkit()`, `FeaturePlot_bioToolkit()`,
  `DensityDiffPlot_bioToolkit()`, `format_sc_plots()`: new `arrow.panels`
  ("bottom.left" (default), "all", "left", "right", "top", "bottom", "none")
  sets which panels get the dimension arrows. In `format_sc_plots()` it
  replaces `add.arrow` and `bottom.left.only` (accepted with a deprecation
  warning), defaults to "bottom.left" (was every plot), and arrows a plot
  already has are replaced rather than added to.
* `VlnPlot_bioToolkit()`: with a lowered `y.max`, violins and boxes are cut
  at the top of the axis (they ran past it) and points above it are not
  drawn.
* `plot_composition_scatter()`: new `error` ("sem" (default), "sd" or
  "ci95", t-based) for the per-sample error bars, now drawn with caps, 0.75 pt
  and black. The source data add `x_sd`, `y_sd`, `x_error`, `y_error`.
* `plot_fc_quadrant()`: the x = 0 and y = 0 lines are 1 pt black; the
  correlation is bold, black and larger; new `trend` adds a least-squares line
  with its 95% confidence band.
* `plot_enrichment_dots()`: strips of all panels have the same thickness,
  names centred; with several panels the legends are in two columns (colour
  bars, then size legends; new `legend.layout`, "stacked" for the previous
  single column); one panel shows its two legends side by side.
* `plot_gsea_rank()`: NES and FDR bold and larger; the gene set ID is the
  second line of the title and the comparison a bold italic subtitle; new
  `label` / `label.n` / `label.features` name the leading-edge (or chosen)
  genes with an arrow to their ticks; new `band = "none"`; new
  `italic.genes`. The source data add `leading_edge` and `labelled`.
* `ggplotColours()` is renamed `ggplot_colors()` (`showColor` ->
  `show.colors`); the old name still works, with a deprecation warning.
* `plot_multi_volcano()`: with `layout = "sina"`, non-significant genes are
  spread uniformly over the column, a grey block centred on 0, instead of
  forming a third shape between the cutoffs; up- and down-regulated genes keep
  their own shapes. With `layout = "pvalue"` (suited to bulk data), each side
  of 0 is scaled as a whole, so grey genes sit at their own significance.
  Non-significant genes are a lighter grey ("grey88"; `plot_volcano()` keeps
  "grey75") and always drawn under the significant ones.
* `plot_gsea_rank()`: leading-edge gene names at the small text size (9 pt,
  8 pt when they need more than two lines; were 11 / 9 pt).
* `show_palettes()` also draws your own colours: a vector of colours, or a
  list of palette names and colour vectors (one row each, labelled by the
  list's names). Its first argument is renamed `palettes` (`names` still
  works, with a deprecation warning).

## Documentation website

* New "Colour palettes" page: every palette with the index and hex code of
  each colour, its type and source, and the plots that use it by default.
* A custom "page not found" page (`.github/404.md`).
* Home page: a general description of the package; plot conventions moved to
  the conventions table.
* Figure essentials: theme tabs fixed and self-contained, panel versus cell
  sizes explained with examples, more legend extraction examples.
* Guides show the new arguments; their code follows the new code style (one
  ggplot layer per line, one argument per line for long calls; see
  CONTRIBUTING.md).

# bioToolkit 0.10.0

Package clean-up ahead of the documentation website
(<https://j-y26.github.io/bioToolkit/>).

## Removed

* Legacy functions no longer maintained: `peak_category_qc()`,
  `peak_category_plot()`, `peak_tss_enrichment()`, `peak_width_distribution()`,
  `peak_calc_frip()`, `peak_plot_frip()`, `grouped_bar_plot()`,
  `calc_module_scores()` and `calc_module_scores_ssgsea()`. Their Suggests
  (ChIPseeker, GenomicRanges, IRanges, Rsamtools, rtracklayer, GSVA,
  BiocParallel, Biobase) and the unused aplot, data.table and stringr were
  dropped. `format_gene_sets()` is kept (now in `R/gene_sets.R`).

## Deprecated argument names

The remaining older functions now use the package's dot.case argument names.
The old names are still accepted, with a deprecation warning (class
`deprecatedWarning`), and positional calls are unchanged:

* `format_sc_plots()`, `manipulate_sc_plot()`, `add_dimension_arrows()`:
  snake_case to dot.case (`remove_axes` -> `remove.axes`, `arrow_length` ->
  `arrow.length`, `legend_position` -> `legend.position`, ...), and `byrow` ->
  `by.row` as in `DimPlot_bioToolkit()`.
* `filter_de_genes()`, `rank_de_genes()`, `generate_de_rank_file()`:
  `de_results` -> `data`, `de_method` -> `method`, `pvalue_col` -> `padj.col`,
  `log2fc_col` -> `logfc.col`, `log2fc_min` -> `logfc.cutoff`, `pvalue_max` ->
  `padj.cutoff`, `save_file` -> `file`, as in `standardize_de()` and
  `plot_volcano()`.

## Changes

* `filter_de_genes()`, `rank_de_genes()`, `generate_de_rank_file()`: `method`
  now defaults to "auto" and also reads limma and Seurat tables (the columns
  come from the presets of `standardize_de()`); the method name ignores case
  as before. `filter_de_genes()` no longer returns all-NA rows for genes with
  an NA adjusted p-value.
* `format_sc_plots()` keeps the relative widths and heights of a patchwork
  (they were passed to `wrap_plots()` under the wrong names and lost).
* `DimPlot_bioToolkit()`, `FeaturePlot_bioToolkit()` and
  `DensityDiffPlot_bioToolkit()` take the dot.case arrow arguments of
  `add_dimension_arrows()` in `...` (e.g. `arrow.length`); the snake_case ones
  are deprecated.
* `format_gene_sets(to = "gmt")` returns the GMT lines (it printed them and
  returned NULL); new `file` argument to write them. A gene listed twice in a
  set (e.g. from msigdbr) is kept once instead of failing.

## Fixes found while writing the website guides

* `plot_multi_volcano()` with `groups`: the comparison names under the
  grouped panels were dropped (out of the y range); they are drawn again.
* `plot_enrichment_dots()`: the legends of stacked panels overlapped the next
  panel's; legends are now collected into one column on the right. Panels
  side by side with different numbers of gene sets had names misaligned with
  their dots; they now have the same height (empty rows at the bottom of
  shorter panels). Long panel names are wrapped to fit their strip.
* `plot_composition_scatter()`: groups with no cells in a panel of
  `facet.by` (or among `cells`) are no longer drawn and labelled at the
  origin.

## Website

* Documentation website with guides by task, run on public data (Kang et
  al. 2018 PBMCs; Alasoo et al. 2018 macrophages; MSigDB Hallmark), and the
  function index in 9 categories. Built by GitHub Actions when the version
  changes. Website files are excluded from the package build, and the
  packages needed only by the guides are in `Config/Needs/website`.

# bioToolkit 0.9.0

Phase E of the visualization overhaul: bulk plots (B1, B2, B4, B5). The
planned sample-to-sample heatmap (B3) was dropped; use corrplot or
ComplexHeatmap.

## New

* `plot_ma()`: MA plot (log2 fold change against mean expression) with the
  classes, colours, labels and arguments of `plot_volcano()`. Mean expression
  as the DE tool reports it: log2 values (edgeR logCPM, limma AveExpr) as
  they are, linear ones (DESeq2 baseMean) on a log10 axis (`mean.scale`);
  optional robust trend line (`trend`, lowess); counts of up / down features
  in the right corners. Shrunken fold changes (DESeq2 `lfcShrink()`) are the
  recommended input (`logfc.col`).
* `plot_pca()`: PCA score plot of samples from a `prcomp` result or a
  normalized matrix (top `n.top` = 500 most variable features), coloured
  (discrete or continuous), shaped and labelled by sample covariates, % of
  variance in the axis titles, optional ellipses (groups of at least 4
  samples) and scree panel; warns on raw counts; square panel; scores and
  covariates as source data.
* `plot_variance_explained()`: variance explained by sample covariates.
  View "genes": per-gene fractions as violins (ordered by median, Residuals
  last), or stacked bars for selected `features`; read from a
  variancePartition result, or computed with linear models (R-squared split
  between covariates by the LMG method, so fractions sum to 1 and do not
  depend on the order of the covariates). View "pcs": R-squared of each top
  PC on each covariate as tiles, under bars of each PC's % variance.
* `plot_fc_quadrant()`: log2 fold change in one comparison against another
  (two time points, two tissues, RNA vs ATAC / ChIP), coloured by where each
  feature is significant (both, same or opposite direction; one only;
  neither), with counts per quadrant, the Spearman (or Pearson) correlation,
  a y = x diagonal and labels of the top features per quadrant. The two
  tables may come from different tools; peaks are mapped to genes with
  `map.y` (one peak per gene: the most significant by default).

## Internal

* `plot_volcano()` uses the shared `check_de_cutoffs()` and the new
  `de_direction()`, also used by `plot_ma()` and `plot_fc_quadrant()`.
  Output unchanged.

# bioToolkit 0.8.0

Phase D of the visualization overhaul: single-cell plots (S1-S4).

## Breaking changes

* `group.by` / `split.by` now mean the same in every plot (user decision
  2026-10-13): `group.by` is the lower-level grouping whose levels are
  compared, marked by their colours; `split.by` is the higher-level grouping,
  drawn like facets (labelled blocks or panels). `DotPlot_bioToolkit()`
  follows it: with `split.by`, rows are blocks of `split.by` (named on the
  left), each with one row per `group.by` level marked by its colour. What was
  `DotPlot_bioToolkit(group.by = "cell_type", split.by = "condition")` is now
  `DotPlot_bioToolkit(group.by = "condition", split.by = "cell_type")` (and
  the dots show the condition colours). This is the opposite of
  `Seurat::DotPlot()` / `VlnPlot()`'s `split.by`, and matches the panels of
  `DimPlot_bioToolkit()` / `FeaturePlot_bioToolkit()`.
* `plot_meta_scatter()` is removed; use `plot_composition_scatter()`.
* `summarize_group_meta()`: `cells.use` is now `cells`; `group.by` is a single
  column; new `facet.by` and `sample.by`; the result has columns `facet`,
  `sample` (when given), `condition`, `group`, `count`, `total` and `prop`
  (was `total_count_condition`). It no longer calls `gc()`.

## New

* `VlnPlot_bioToolkit()`: per-cell violins with a box (quartiles and median)
  and a random sample of cells as rasterized points, for genes, metadata (QC
  metrics, module scores) or both, from a Seurat object or a data frame. One
  panel (and value axis) per feature: `log` and `percent` take TRUE or feature
  names; `y.max` zooms to a number or a quantile ("q99") without dropping
  cells; `zeros = "exclude"` / `"annotate"` for sparse genes (violins of the
  cells above 0, with the % above 0); `flip` (horizontal), `stack` (thin rows
  sharing the group axis); `split.by` as labelled blocks of violins; fixed
  violin sizes; source data per cell with the assay and layer.
* `RidgePlot_bioToolkit()`: density ridges with the same data arguments
  (needs ggridges, now in Suggests); `overlap` sets the ridge height; one
  bandwidth per feature.
* `DoHeatmap_bioToolkit()`: heatmap of genes x cells (`columns = "cells"`,
  at most `max.cells` per group, seeded) or group means (`columns =
  "groups"`), z-scored per gene from the "data" layer by default (across all
  plotted cells, before downsampling, or across the columns), so that any
  gene can be drawn (Seurat's "scale.data" holds only variable genes).
  Columns in blocks of `split.by` (named, with gaps) and groups compared
  within them, marked by a coloured dot (default) or strip with a legend;
  gene blocks from a named list (e.g. `top_markers()`); extra annotation
  strips (`annotations`); gene names drawn when rows are tall enough, on the
  right or left (`label.side`), or a few named with leader lines
  (`label.features`). No clustering (use ComplexHeatmap).
* `plot_composition_scatter()`: proportion (or count) of each group in one
  condition against another, with the y = x line: `split.by` with levels `x`
  and `y` on the axes, `facet.by` panels (proportions within each), optional
  per-sample proportions with mean +- SEM (`sample.by`; pooled by default),
  `log` axes, `fold.lines`, labels for all or selected groups; square panels
  with the same range on both axes.
* `DensityDiffPlot_bioToolkit()`: where on an embedding the cells of a
  condition are enriched or depleted relative to a reference. One kernel
  density per condition (one bandwidth for all), normalized to sum to 1 so
  that cell numbers do not matter; optionally one density per sample,
  averaged per condition (`sample.by`); `measure = "difference"` (relative to
  the reference's peak) or `"log2ratio"`; sparse regions masked; cells
  coloured by the value at their position (`draw = "cells"`) or the grid as a
  raster; one panel per compared condition, one symmetric diverging scale;
  optional density contours, group labels and dimension arrows.
* Palette "teal-yellow" (teal - black - yellow, built as Seurat's
  PurpleAndYellow), the default of `DoHeatmap_bioToolkit()` for z-scores of
  cells; heatmaps of group means default to "blue-magenta" (light centre).
* Genomic regions (e.g. peaks "chr1-100-200") are no longer set in italics as
  gene names.
* Internal: the kernel smoothing of `FeaturePlot_bioToolkit(geom = "density")`
  is split into `kde_grid()` and `grid_interpolate()`, shared with
  `DensityDiffPlot_bioToolkit()` (results unchanged).

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
