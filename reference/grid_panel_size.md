# Panel size of a grid of categories

With \`cell\` inches per category, a panel holding \`n\` categories is
\`n \* cell\` inches. For facets of different numbers of categories
(e.g. \`facet_grid(space = "free")\`), give \`n\` per facet column (or
row) to get one panel size per column (row), for
[`panel_size_theme`](https://j-y26.github.io/bioToolkit/reference/panel_size_theme.md).

## Usage

``` r
grid_panel_size(n, cell)
```

## Arguments

- n:

  Number of categories, one value per panel column (or row)

- cell:

  Size of one category in inches

## Value

Panel sizes in inches, one per element of \`n\`
