# Restyle a single ggplot made by Seurat

Removes the axes, legend or grid of one ggplot (e.g. from
[`Seurat::DimPlot()`](https://satijalab.org/seurat/reference/DimPlot.html)),
restores a hidden legend, or sets the legend position and title. Used by
[`format_sc_plots`](https://j-y26.github.io/bioToolkit/reference/format_sc_plots.md)
on each plot of a list; only valid for individual ggplot objects.

## Usage

``` r
manipulate_sc_plot(
  plot,
  remove.axes = FALSE,
  remove.legend = FALSE,
  remove.grid = FALSE,
  restore.legend = FALSE,
  legend.position = NULL,
  legend.title = NULL,
  legend.title.size = 12,
  legend.title.bold = TRUE,
  ...
)
```

## Arguments

- plot:

  A single ggplot object

- remove.axes:

  Whether to remove the axes, default is FALSE

- remove.legend:

  Whether to remove the legend, default is FALSE

- remove.grid:

  Whether to remove the grid, default is FALSE

- restore.legend:

  Whether to restore the legend, default is FALSE

- legend.position:

  The position of the legend, default is NULL (meaning no change)

- legend.title:

  The title of the legend, default is NULL

- legend.title.size:

  The size of the legend title, default is 12

- legend.title.bold:

  Whether to bold the legend title, default is TRUE

- ...:

  Deprecated snake_case names of the arguments above (e.g.
  \`remove_axes\`), accepted with a warning; renamed to dot.case in
  0.10.0

## Value

a ggplot object with the axes, legend, and grid manipulated
