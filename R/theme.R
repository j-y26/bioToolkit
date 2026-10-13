# Shared theme and text helpers (roadmap item F3)
#
# Design
# ------
# * One theme, `theme_bioToolkit()`, for every plot in the package, so that
#   text sizes, line widths and fonts are consistent across figures.
# * One fixed set of sizes, `theme_sizes()`, used everywhere. (Journal and
#   slide presets existed briefly and were removed on 2026-10-13 to keep the
#   package simple.) Each plot function still accepts its own `*.size`
#   arguments, whose NULL default means "from `theme_sizes()`".
# * The sizes, all in points (1/72 inch):
#     base  - axis titles, legend titles and, in embeddings, the legend text
#             (cluster names are the content of an embedding legend, so they
#             are not demoted to "small" text)
#     small - tick labels and legend text of axis / grid plots
#     large - plot titles and facet strips
#     line  - thin lines: reference lines (e.g. cutoffs), label leader
#             segments and grid lines
#     axis  - axis lines and ticks of "axis" and "grid" plots, firmer than the lines
#             drawn inside the panel (user decision 2026-10-13)
# * ggplot2 measures `linewidth` in its own units, so line widths given in
#   points are converted with `pt_to_linewidth()`; `geom_text()` sizes are in
#   mm and converted with `pt_to_size()`.


#' Text sizes and line widths of `theme_bioToolkit()`
#'
#' @param base.size If not NULL, the base text size; the small and large sizes
#'                  are scaled with it so that their ratios are kept
#'
#' @return A list with `base` (11), `small` (9), `large` (12), `line` (0.5;
#'         thin lines) and `axis` (1; axis lines and ticks), in points
#'
#' @keywords internal
theme_sizes <- function(base.size = NULL) {
    sizes <- list(base = 11, small = 9, large = 12, line = 0.5, axis = 1)
    if (!is.null(base.size)) {
        ratio <- base.size / sizes$base
        sizes$small <- sizes$small * ratio
        sizes$large <- sizes$large * ratio
        sizes$base <- base.size
    }
    sizes
}


#' Convert a line width in points to ggplot2 `linewidth` units
#'
#' ggplot2 draws a `linewidth` of 1 as `ggplot2::.pt` device line widths, each
#' of which is 1/96 inch, i.e. 0.75 pt.
#'
#' @param pt Line width in points
#'
#' @return The `linewidth` giving that width
#'
#' @keywords internal
pt_to_linewidth <- function(pt) {
    pt * 96 / 72 / ggplot2::.pt
}


#' Convert a line width in points to the `stroke` of ggplot2 points
#'
#' ggplot2 draws the outline of a point `stroke * ggplot2::.stroke / 2`
#' device line widths wide, each 1/96 inch, i.e. 0.75 pt.
#'
#' @param pt Outline width in points
#'
#' @return The `stroke` giving that width
#'
#' @keywords internal
pt_to_stroke <- function(pt) {
    pt * 96 / 72 * 2 / ggplot2::.stroke
}


#' Convert a font size in points to the mm used by `geom_text()` `size`
#'
#' @param pt Font size in points
#'
#' @keywords internal
pt_to_size <- function(pt) {
    pt / ggplot2::.pt
}


#' @title The bioToolkit ggplot2 theme
#'
#' @description The theme shared by all plots of the package: 11 pt text
#'              (9 pt tick labels and legend text, 12 pt titles), 1 pt axis
#'              lines and 0.5 pt other lines. Three layouts are provided:
#'              \describe{
#'                \item{"axis"}{Classic axes: axis lines and black ticks, no
#'                      grid. For scatter, volcano, MA, bar and PCA plots.}
#'                \item{"grid"}{The axes of "axis", plus light grid lines,
#'                      no axis titles and x labels rotated by 45 degrees. For
#'                      dot plots and heatmaps, whose rows and columns are
#'                      categories.}
#'                \item{"embedding"}{No axes at all, with extra room at the
#'                      bottom left for the dimension arrows. For UMAP / tSNE
#'                      plots.}
#'              }
#'
#'              Panel sizes are not part of the theme; add
#'              \code{\link{panel_size_theme}} to fix them.
#'
#' @param type Layout of the theme, one of "axis" (default), "grid" or
#'             "embedding"
#'
#' @param base.size Base font size in points. If NULL (default), 11; otherwise
#'                  the other text sizes are scaled with it
#'
#' @param base.family Font family. If NULL (default), the device default
#'                    (Helvetica in PDF, which journals accept as Arial's
#'                    equivalent)
#'
#' @param bg Background colour, default is "transparent" (NA also means
#'           transparent), so that figures can be placed on any background.
#'           Use e.g. "white" for an opaque background
#'
#' @return A ggplot2 theme
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
#' p + theme_bioToolkit()
#' p + theme_bioToolkit(base.size = 8)
#'
#' @seealso \code{\link{panel_size_theme}}, \code{\link{save_plot}},
#'          \code{\link{italic_labels}}
#'
#' @import ggplot2
#' @export
#'
theme_bioToolkit <- function(type = c("axis", "grid", "embedding"),
                             base.size = NULL,
                             base.family = NULL,
                             bg = "transparent") {
    type <- match.arg(type)
    sz <- theme_sizes(base.size)
    family <- if (is.null(base.family)) "" else base.family
    lw <- pt_to_linewidth(sz$line)
    transparent_bg <- is.na(bg) || identical(bg, "transparent")
    plot_bg <- element_rect(fill = if (transparent_bg) "transparent" else bg, colour = NA)

    # Elements shared by the three layouts: legends without a box, strips
    # without a background, and the plot background
    common <- theme(
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        # Gap between panel and legend on whichever side the legend is
        legend.box.spacing = unit(0.3, "cm"),
        strip.background = element_blank(),
        plot.background = plot_bg,
        panel.background = element_blank(),
        plot.title = element_text(size = sz$large, face = "bold", hjust = 0,
                                  margin = margin(b = 0.15, unit = "cm")),
        plot.title.position = "panel"
    )

    if (type == "embedding") {
        # Margins scale with a custom base size (1 at the default 11 pt)
        k <- sz$base / theme_sizes()$base
        # The look of the reduction plots: nothing but the cells, the legend
        # and the strips. The legend text is the base size because it names the
        # groups, which are the content of the plot
        return(
            theme_void(base_size = sz$base, base_family = family) +
                common +
                theme(
                    legend.text = element_text(size = sz$base),
                    legend.title = element_text(size = sz$base),
                    strip.text = element_text(size = sz$large, face = "bold",
                                              margin = margin(b = 0.15, unit = "cm")),
                    # Extra bottom / left margin makes room for the arrow labels;
                    # it scales with a custom base size
                    plot.margin = margin(0.2 * k, 0.2 * k, 0.6 * k, 0.6 * k, "cm")
                )
        )
    }

    th <- theme_classic(base_size = sz$base, base_family = family,
                        base_line_size = lw, base_rect_size = lw) +
        common +
        theme(
            axis.title = element_text(size = sz$base, colour = "black"),
            axis.text = element_text(size = sz$small, colour = "black"),
            # Axes use their own, firmer width (`axis`)
            axis.line = element_line(colour = "black", linewidth = pt_to_linewidth(sz$axis)),
            axis.ticks = element_line(colour = "black", linewidth = pt_to_linewidth(sz$axis)),
            # Ticks scale with the text, so they stay in proportion with a
            # custom base size
            axis.ticks.length = unit(sz$small / 3, "pt"),
            legend.title = element_text(size = sz$base),
            legend.text = element_text(size = sz$small),
            strip.text = element_text(size = sz$base, face = "bold",
                                      margin = margin(b = 0.1, t = 0.1, unit = "cm")),
            plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")
        )

    if (type == "grid") {
        # Rows and columns are categories: the 1 pt axes of "axis" plots
        # (user decision 2026-10-13; they were dropped before), light grid
        # lines to follow a row or column, no axis titles (the labels name
        # the categories) and rotated x labels
        th <- th + theme(
            axis.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            panel.grid.major = element_line(colour = "grey92", linewidth = lw)
        )
    }
    th
}


#' @title Labels with gene names in italics
#'
#' @description Gene symbols are set in italics by convention. This turns
#'              labels into plotmath expressions, italic where they are genes
#'              and plain otherwise, and can be passed directly as the `labels`
#'              of a discrete scale. Labels are used as they are: the case of
#'              gene symbols (e.g. "CD3E" in human, "Cd3e" in mouse and rat) is
#'              not changed.
#'
#' @param x Character vector of labels
#'
#' @param is.gene Whether each label is a gene: a single TRUE / FALSE for all,
#'                a logical vector along `x`, or a logical vector named by
#'                label. Labels missing from a named vector are plain
#'
#' @return An expression vector, one element per label
#'
#' @examples
#' library(ggplot2)
#' df <- data.frame(gene = c("Cd3e", "Ms4a1", "nCount_RNA"), value = 1:3)
#' ggplot(df, aes(value, gene)) +
#'     geom_point() +
#'     scale_y_discrete(labels = function(x) italic_labels(x, x != "nCount_RNA"))
#'
#' @export
#'
italic_labels <- function(x,
                          is.gene = TRUE) {
    x <- as.character(x)
    if (!is.null(names(is.gene))) {
        flag <- is.gene[x]
        is.gene <- !is.na(flag) & flag
    } else {
        is.gene <- rep_len(as.logical(is.gene), length(x))
    }
    # bquote() keeps any character in a label, where parse() would fail on
    # symbols such as "-" or "'"
    as.expression(lapply(seq_along(x), function(i) {
        if (isTRUE(is.gene[i])) bquote(italic(.(x[i]))) else x[i]
    }))
}
