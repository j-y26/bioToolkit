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
