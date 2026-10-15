# Colours of a palette

Returns the colours of a palette of the package as a plain vector of hex
codes, to pick from, reorder, name and pass to the \`cols\` argument of
any plot function (or to `scale_fill_manual()`, ComplexHeatmap, base
graphics...). Use
[`show_palettes`](https://j-y26.github.io/bioToolkit/reference/show_palettes.md)
to see the colours with their index.

## Usage

``` r
bt_palette(name, n = NULL, reverse = FALSE)
```

## Arguments

- name:

  Palette name (case is ignored), see
  [`list_palettes`](https://j-y26.github.io/bioToolkit/reference/list_palettes.md).
  Add "\_r" to any name for the palette reversed, e.g. "viridis_r"
  (yellow to dark purple) or "npg_r"

- n:

  Number of colours. If NULL (default), all colours of the palette.
  Qualitative palettes return their first \`n\` colours and are never
  interpolated, since a colour in between two group colours would read
  as an intermediate group; sequential and diverging palettes are
  interpolated to \`n\` colours

- reverse:

  Whether to reverse the order of the colours, default is FALSE

## Value

A character vector of hex colours

## See also

[`show_palettes`](https://j-y26.github.io/bioToolkit/reference/show_palettes.md),
[`list_palettes`](https://j-y26.github.io/bioToolkit/reference/list_palettes.md)

## Examples

``` r
bt_palette("scanpy-20")
#>  [1] "#1F77B4" "#FF7F0E" "#2CA02C" "#D62728" "#9467BD" "#8C564B" "#E377C2"
#>  [8] "#7F7F7F" "#BCBD22" "#17BECF" "#AEC7E8" "#FFBB78" "#98DF8A" "#FF9896"
#> [15] "#C5B0D5" "#C49C94" "#F7B6D2" "#C7C7C7" "#DBDB8D" "#9EDAE5"

# Pick and reorder colours by index, then name them by cell type
cols <- bt_palette("scanpy-20")[c(1, 3, 2, 10)]
names(cols) <- c("T cells", "B cells", "Monocytes", "NK cells")
cols
#>   T cells   B cells Monocytes  NK cells 
#> "#1F77B4" "#2CA02C" "#FF7F0E" "#17BECF" 

# Low, mid and high colours of a diverging scale
bt_palette("tol-burd")[c(1, 5, 9)]
#> [1] "#2166AC" "#F7F7F7" "#B2182B"
```
