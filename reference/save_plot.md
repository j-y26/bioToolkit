# Save a plot to a file sized to fit it

Saves a ggplot or patchwork to PDF (vector), PNG or any format
[`ggsave`](https://ggplot2.tidyverse.org/reference/ggsave.html)
supports.

For plots whose panels have a fixed size (all plots of this package, or
any ggplot with
[`panel_size_theme`](https://j-y26.github.io/bioToolkit/reference/panel_size_theme.md)),
the width and height of the file are computed from the panels, strips,
legends and margins, so that nothing is cropped and no white space is
added. They are measured on the same kind of device the file is written
with, so that font metrics match. Patchworks of such plots, including
nested ones, keep the size of every panel.

Plots without fixed panel sizes are saved at \`width\` by \`height\`
inches, as with
[`ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

A legend from
[`extract_legend`](https://j-y26.github.io/bioToolkit/reference/extract_legend.md)
is saved on its own, sized to fit, with a small margin.

Files are written with cairo when R has it (`capabilities("cairo")`):
[`cairo_pdf`](https://rdrr.io/r/grDevices/cairo.html) for PDFs, and
[`png`](https://rdrr.io/r/grDevices/png.html) (or
[`jpeg()`](https://rdrr.io/r/grDevices/png.html),
[`tiff()`](https://rdrr.io/r/grDevices/png.html),
[`bmp()`](https://rdrr.io/r/grDevices/png.html)) with `type = "cairo"`
for raster files. Fonts are embedded, Unicode text is drawn, a hyphen
stays a hyphen (base [`pdf()`](https://rdrr.io/r/grDevices/pdf.html)
draws "-" with the wider minus glyph), and PDF and PNG share the same
font. Without cairo, for other formats (e.g. SVG), or when \`device\` is
given,
[`ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)'s
device is used and \`family\` is ignored.

## Usage

``` r
save_plot(
  plot,
  filename,
  width = NULL,
  height = NULL,
  dpi = 1000,
  bg = NULL,
  family = "Arial",
  ...
)
```

## Arguments

- plot:

  A ggplot or patchwork object, or a legend from
  [`extract_legend`](https://j-y26.github.io/bioToolkit/reference/extract_legend.md)

- filename:

  File to save to; the format is inferred from the extension (".pdf",
  ".png", ".svg", ".tiff", ...)

- width, height:

  Size of the file in inches, used only for plots without fixed panel
  sizes, for which they are required. Default is NULL

- dpi:

  Resolution of raster formats, default is 1000

- bg:

  Background colour. If NULL (default), taken from the plot. A
  transparent background is kept in PDF, PNG, TIFF and SVG; formats
  without transparency (JPEG, BMP) get a white one instead

- family:

  Font family of all text that has no font of its own (the default of
  ggplot2 and of
  [`theme_bioToolkit()`](https://j-y26.github.io/bioToolkit/reference/theme_bioToolkit.md);
  a \`base.family\` set in the theme is kept). Default is "Arial". If
  the font is not installed, the system substitutes one (e.g. Liberation
  Sans, which has Arial's widths, or the wider DejaVu Sans). NULL keeps
  the devices' default "sans"

- ...:

  Additional arguments passed to
  [`ggsave`](https://ggplot2.tidyverse.org/reference/ggsave.html)

## Value

The width and height of the saved file in inches, invisibly

## See also

[`panel_size_theme`](https://j-y26.github.io/bioToolkit/reference/panel_size_theme.md),
[`export_source_data`](https://j-y26.github.io/bioToolkit/reference/export_source_data.md),
[`extract_legend`](https://j-y26.github.io/bioToolkit/reference/extract_legend.md)

## Examples

``` r
library(ggplot2)
p <- ggplot(mtcars, aes(wt, mpg)) + geom_point() +
    theme_bioToolkit() + panel_size_theme(2.5, 2)
f <- tempfile(fileext = ".pdf")
save_plot(p, f)
```
