# Shared setup of the bioToolkit website guides, sourced by each guide's
# hidden setup chunk. Files starting with "_" are not articles for pkgdown.

suppressPackageStartupMessages({
    library(bioToolkit)
    library(ggplot2)
})
source("_data.R")

knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    # Package messages (e.g. which assay / DE method was used) are shown only
    # in the chunks that explain them
    message = FALSE,
    fig.align = "center"
)

# Figures ------------------------------------------------------------------
# Every ggplot, patchwork or extracted legend printed in a guide is written
# with save_plot(), at the size the plot fits (bioToolkit plots have fixed
# panel sizes), so the web figure is exactly what save_plot() gives a user.
# 192 dpi is sharp on high-density screens; the image is shown at its
# physical size (include_graphics(dpi = 192): 1 inch = 96 CSS pixels).
# Plots without fixed panels use the chunk's fig.width / fig.height.
WEB_DPI <- 192
fig_counter <- new.env()
fig_counter$n <- 0

bt_knit_print <- function(x,
                          options,
                          ...) {
    fig_counter$n <- fig_counter$n + 1
    file <- knitr::fig_path(".png", options, number = fig_counter$n)
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    suppressMessages(bioToolkit::save_plot(x, file, dpi = WEB_DPI,
                                           width = options$fig.width, height = options$fig.height))
    knitr::knit_print(knitr::include_graphics(file, dpi = WEB_DPI), options)
}
for (cls in c("ggplot", "patchwork", "bt_legend")) {
    registerS3method("knit_print", cls, bt_knit_print, envir = asNamespace("knitr"))
}
