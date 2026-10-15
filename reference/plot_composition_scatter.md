# Composition scatter: group proportions in two conditions

Plots, for each group of cells (e.g. cell type or state), its proportion
(or number of cells) in one condition on x against another condition on
y, with the y = x line: groups above the line are more abundant in the y
condition. Optionally one panel per level of \`facet.by\` (e.g. stage),
and per-sample proportions with their mean and SEM (\`sample.by\`).
Square panels of a fixed size; save it with
[`save_plot`](https://j-y26.github.io/bioToolkit/reference/save_plot.md);
the table behind it is returned by
[`export_source_data`](https://j-y26.github.io/bioToolkit/reference/export_source_data.md).

## Usage

``` r
plot_composition_scatter(
  object,
  group.by,
  split.by,
  x = NULL,
  y = NULL,
  facet.by = NULL,
  sample.by = NULL,
  error = c("sem", "sd", "ci95"),
  cells = NULL,
  groups = NULL,
  value = c("percent", "count"),
  log = FALSE,
  fold.lines = NULL,
  label = TRUE,
  label.size = NULL,
  cols = NULL,
  pt.size = 4,
  panel.width = 2.5,
  panel.height = panel.width,
  legend.position = "none",
  bg = "transparent"
)
```

## Arguments

- object:

  A Seurat object, or a data frame with one row per cell (e.g. its
  metadata)

- group.by:

  Column of the groups, one point each (e.g. cell state)

- split.by:

  Column of the conditions; its levels \`x\` and \`y\` are compared on
  the axes (e.g. genotype)

- x, y:

  Levels of \`split.by\` on the x- and y-axis. If NULL (default), the
  first and second level

- facet.by:

  Column giving one panel per level (e.g. stage), with proportions
  computed within each. Default is NULL

- sample.by:

  Column of the samples. If given, proportions are computed per sample
  and each point is the mean of the samples of a condition, with error
  bars on both axes (see \`error\`), so that every sample counts the
  same. If NULL (default), cells of all samples are pooled into one
  proportion per condition

- error:

  Error bars with \`sample.by\`: "sem" (default; standard error of the
  mean), "sd" (standard deviation of the samples) or "ci95" (95%
  confidence interval of the mean, mean +- t(0.975, n - 1) x SEM). The
  horizontal bar is computed from the samples of the x condition, the
  vertical one from those of the y condition. A condition with one
  sample has no bar

- cells:

  Cells (rows) to use, as names or indices, e.g. one lineage:
  proportions are then of that subset. If NULL (default), all

- groups:

  Groups to draw. Proportions are still computed over all groups. If
  NULL (default), all

- value:

  "percent" (default; the proportion, labelled in %) or "count" (the
  number of cells)

- log:

  Whether to draw both axes on a log10 scale, default is FALSE. Spreads
  rare groups; groups with 0 cells in a condition cannot be drawn and
  are left out, with a message

- fold.lines:

  Fold changes drawn as dashed lines around y = x, e.g. 2 draws y = 2x
  and y = x / 2. Default is NULL (none)

- label:

  Which points to name (needs ggrepel): TRUE (default) for all, FALSE
  for none, or the groups to name. Labels of many rare groups crowd near
  the origin on linear axes; name a few, or use \`log = TRUE\`

- label.size:

  Size of the labels in mm. If NULL (default), the small text size of
  the theme (9 pt)

- cols:

  Colours of the groups: NULL (default) for the colours
  [`DimPlot_bioToolkit`](https://j-y26.github.io/bioToolkit/reference/DimPlot_bioToolkit.md)
  gives them, a qualitative palette name, or colours named by group

- pt.size:

  Size of the points, default is 4

- panel.width, panel.height:

  Size of each (square) panel in inches, default 2.5; \`panel.height\`
  defaults to \`panel.width\`

- legend.position:

  Position of the group legend, default is "none" (points are named)

- bg:

  Background colour, default is "transparent"

## Value

A ggplot object

## Details

The source data
([`export_source_data`](https://j-y26.github.io/bioToolkit/reference/export_source_data.md))
hold one row per group (and panel): the values on x and y, the cell
counts and totals of both conditions and, with \`sample.by\`, the SEM,
SD, number of samples and the half-width of the error bars drawn
(\`x_error\`, \`y_error\`). Per-sample proportions are returned by
[`summarize_group_meta`](https://j-y26.github.io/bioToolkit/reference/summarize_group_meta.md).

Groups with no cells in either condition (of a panel, or among
\`cells\`) are not drawn.

The plot shows where proportions differ, not whether the difference is
significant: that needs a test on samples, not cells (e.g. propeller,
scCODA).

## See also

[`summarize_group_meta`](https://j-y26.github.io/bioToolkit/reference/summarize_group_meta.md),
[`save_plot`](https://j-y26.github.io/bioToolkit/reference/save_plot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Pooled proportions of each cell state, WT vs Mut, one panel per stage
plot_composition_scatter(obj, group.by = "cell_state", split.by = "genotype",
                         x = "WT", y = "Mut", facet.by = "stage")

# Per-sample proportions: mean +- SEM; log axes for rare states
p <- plot_composition_scatter(obj, group.by = "cell_state", split.by = "genotype",
                              x = "WT", y = "Mut", sample.by = "sample",
                              log = TRUE, fold.lines = 2)
save_plot(p, "composition.pdf")
} # }
```
