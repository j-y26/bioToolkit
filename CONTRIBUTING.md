# Contributing to bioToolkit

## Naming conventions

One naming scheme across the package. Every new function follows it; existing
functions are renamed when their roadmap item is worked on.

### Function names

| Kind | Scheme | Examples |
|---|---|---|
| Plots of a Seurat object, per cell | `<SeuratName>_bioToolkit` (drop-in replacements users recognise) | `DimPlot_bioToolkit`, `FeaturePlot_bioToolkit`, `DotPlot_bioToolkit`, `VlnPlot_bioToolkit`, `DoHeatmap_bioToolkit`, `DensityDiffPlot_bioToolkit` |
| Plots of tables / matrices (DE, enrichment, bulk) | `plot_<thing>()` | `plot_volcano`, `plot_ma`, `plot_pca`, `plot_sample_heatmap`, `plot_variance_explained`, `plot_fc_quadrant`, `plot_de_dots`, `plot_enrichment_dots`, `plot_enrichment_bars`, `plot_composition_scatter` |
| Data helpers | verb first | `standardize_de()`, `standardize_enrichment()`, `top_markers()`, `export_source_data()` |
| Palettes | `bt_palette()`; verb first for the rest | `bt_palette()`, `list_palettes()`, `show_palettes()`, `check_palette()` |
| ggplot2 extensions | ggplot2's own scheme, with `bioToolkit` | `theme_bioToolkit()` |
| Saving | `save_plot()` | replaces `save_dimplot()`; `save_umap()` removed |

### Argument names

Seurat-style `dot.case` everywhere, with one shared vocabulary: the same
meaning has the same name and the same default in every function.

| Argument | Meaning |
|---|---|
| `group.by`, `split.by` | grouping / panel variables |
| `cols` | discrete colours (named or ordered) or a palette name — never `colors` / `palette` |
| `na.value` | colour of missing / non-coloured marks |
| `pt.size`, `alpha` | point size and opacity |
| `label`, `label.size`, `label.n` | text labels and how many to draw |
| `legend.position`, `legend.text.size`, `legend.title` | legend |
| `panel.width`, `panel.height` (inches) | fixed panel size |
| `cell.width`, `cell.height` (inches) | fixed size of one row / column of a grid plot |
| `raster`, `raster.dpi` | rasterize the marks only; `raster.dpi` defaults to 1000 everywhere |
| `bg` | background |
| `seed` | randomness, without touching the global RNG (`with_seed()`) |
| `logfc.col`, `p.col`, `padj.col`, `mean.col`, `feature.col` | columns of a DE table |

Default resolution: every rasterized output (`raster.dpi`, the `dpi` of
`save_plot()` and of functions that save figures) defaults to 1000 dpi.

There are no size presets: one set of sizes (`theme_sizes()`: 11 / 9 / 12 pt
text, 1 pt axes, 0.5 pt other lines) applies everywhere. Text-size arguments
(`legend.text.size`, `strip.text.size`, ...) default to `NULL`, meaning "from
`theme_sizes()`", so that the text sizes are set in one place; other sizes
have plain numeric defaults.

## Code style

Function definitions list one parameter per line, never several on a line:

```r
plot_volcano <- function(
    data,
    method = c("auto", "edgeR", "DESeq2", "limma", "seurat"),
    feature.col = NULL
) {
```

or, aligned with the parenthesis:

```r
standardize_de <- function(data,
                           method = c("auto", "edgeR", "DESeq2", "limma", "seurat"),
                           feature.col = NULL) {
```

This applies to all package-level functions, exported or internal. Short
anonymous functions inside code (`function(x) ...`) are exempt.

## Every plot

Before a plotting item is marked done in `dev/visualization_roadmap.md`:

- Follows the names above.
- Fixed geometry: `panel_size_theme()` (continuous axes: `panel.width`,
  `panel.height`, with `panel.height = panel.width` as the default for
  scatter-type plots usually drawn square (volcano, PCA, embeddings, quadrant
  and composition scatter), so that `panel.width` alone resizes them; discrete axes: `cell.width`, `cell.height` with
  `grid_panel_size()` and `grid_expansion()`), and saves correctly with
  `save_plot()`.
- Uses `theme_bioToolkit()`. Colours: one `cols` argument with a good
  default, resolved by `group_colours()` for groups (palette name, or colours
  named by group or in level order) and `resolve_continuous_palette()` for
  continuous scales (palette name, or colours from low to high, e.g.
  `c(low, mid, high)`).
- Gene names in italics (`italic_labels()`, `feature_is_gene()`).
- Expression from Seurat objects: an `assay` argument (NULL = the default
  assay, as in Seurat) and `slot`; resolve them with `resolve_assay()`, which
  reports the assay and layer used, and fetch genes with `fetch_features()`
  (assay-keyed, so values never come from another assay). Write the assay and
  layer into the source data (`R/seurat_data.R`).
- Grids of dots (genes x groups, terms x comparisons) use the dot-grid engine
  `dot_grid_plot()` (`R/dot_grid.R`) rather than drawing their own.
- Legend extractable with `extract_legend()` (works for any ggplot; nothing
  to do unless the plot hides its legend in an unusual way).
- Source data recoverable with `export_source_data()`: the tidy table in
  `plot$data`, or in `attr(plot, "source_data")` when layers use different
  data.
- roxygen docs with an example; `devtools::document()` re-run.
- Unit tests for the data logic in `tests/testthat/`.
- Checked by eye on real data.
- `NEWS.md` entry, including anything removed or renamed.

## Adding a colour palette

Add an entry to `palette_list()` in `R/palettes.R` with its source, copying
hex values from the primary source (or a package reproducing it), never from
memory. The tests in `tests/testthat/test-palettes.R` check every entry.
