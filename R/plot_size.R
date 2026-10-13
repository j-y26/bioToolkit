# Fixed-size plotting engine (roadmap item F2)
#
# Design
# ------
# Every plot of the package has an exact, absolute geometry: its panels are a
# set number of inches wide and tall, whatever the number of panels, the
# legend or the device. Two consequences follow:
#   * panels line up across figures and inside patchworks, and
#   * the file a plot is saved to is computed to fit the plot (`save_plot()`),
#     rather than the plot being stretched to fit a file size chosen by hand.
#
# Panel sizes come from one of two sources:
#   * continuous axes (embeddings, volcano, MA, PCA): `panel.width` and
#     `panel.height` in inches;
#   * discrete axes (dot plots, heatmaps, violins): `cell.width` and
#     `cell.height`, the size of one category. The panel is then
#     `n * cell` inches, so that one gene row is the same height in every
#     figure, whatever the number of genes (`grid_panel_size()`). The scale
#     must then be expanded by exactly half a category on each side
#     (`grid_expansion()`), so that n categories tile the panel exactly.
#
# Both end in `panel_size_theme()`, which only sets `panel.widths`,
# `panel.heights` and `panel.spacing` (ggplot2 >= 4.0.0). The look of a plot
# is set separately by `theme_bioToolkit()` (F3).


#' @title Theme giving panels a fixed absolute size
#'
#' @description Sets the size of every panel of a plot in inches, so that it is
#'              identical across figures whatever the number of panels and the
#'              size or position of the legend. Add it to any ggplot, then save
#'              the plot with \code{\link{save_plot}}, which sizes the file to
#'              fit. Requires ggplot2 >= 4.0.0.
#'
#' @param panel.width Width of a panel in inches: one value for all panel
#'                    columns, or one value per column (e.g. from facets of
#'                    different numbers of categories)
#'
#' @param panel.height Height of a panel in inches: one value for all panel
#'                     rows, or one per row. Default is `panel.width` when it
#'                     is a single value
#'
#' @param ncol,nrow Number of panel columns and rows, default is 1. Ignored
#'                  when `panel.width` (`panel.height`) gives one value per
#'                  column (row)
#'
#' @param panel.spacing Gap between panels in inches. If NULL (default), the
#'                      theme's own spacing is kept
#'
#' @return A ggplot2 theme
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point() +
#'     theme_bioToolkit() + panel_size_theme(2.5, 2)
#' \dontrun{
#' save_plot(p, "mtcars.pdf")
#' }
#'
#' @seealso \code{\link{save_plot}}, \code{\link{theme_bioToolkit}}
#'
#' @export
#'
panel_size_theme <- function(panel.width, panel.height = panel.width, ncol = 1, nrow = 1,
                             panel.spacing = NULL) {
    check_fixed_panel_support()
    check_size <- function(x, arg) {
        if (!is.numeric(x) || length(x) == 0 || anyNA(x) || any(x <= 0)) {
            stop("`", arg, "` must be positive numbers (inches)", call. = FALSE)
        }
    }
    check_size(panel.width, "panel.width")
    check_size(panel.height, "panel.height")
    # Given one unit per panel column / row, ggplot2 sizes each panel rather
    # than the panel area as a whole
    widths <- if (length(panel.width) > 1) panel.width else rep(panel.width, ncol)
    heights <- if (length(panel.height) > 1) panel.height else rep(panel.height, nrow)
    th <- theme(
        panel.widths = unit(widths, "in"),
        panel.heights = unit(heights, "in")
    )
    if (!is.null(panel.spacing)) th <- th + theme(panel.spacing = unit(panel.spacing, "in"))
    th
}


#' Panel size of a grid of categories
#'
#' With `cell` inches per category, a panel holding `n` categories is
#' `n * cell` inches. For facets of different numbers of categories (e.g.
#' `facet_grid(space = "free")`), give `n` per facet column (or row) to get
#' one panel size per column (row), for \code{\link{panel_size_theme}}.
#'
#' @param n Number of categories, one value per panel column (or row)
#'
#' @param cell Size of one category in inches
#'
#' @return Panel sizes in inches, one per element of `n`
#'
#' @keywords internal
grid_panel_size <- function(n, cell) {
    if (!is.numeric(cell) || length(cell) != 1 || is.na(cell) || cell <= 0) {
        stop("Cell sizes must be a single positive number (inches)", call. = FALSE)
    }
    # A facet with no category would get a zero-size panel, which ggplot2
    # rejects; it gets one cell instead
    pmax(n, 1) * cell
}


#' Scale expansion that makes categories tile a grid panel exactly
#'
#' A discrete scale of n categories spans 1 to n; expanding it by half a
#' category on each side makes it span n units, so that each category gets
#' exactly the `cell` size given to \code{\link{grid_panel_size}}. (ggplot2's
#' default of 0.6 would make cells slightly smaller than asked.)
#'
#' @return An expansion vector for the `expand` argument of a discrete scale
#'
#' @keywords internal
grid_expansion <- function() {
    expansion(add = 0.5)
}


#' Whether a gtable has an absolute size
#'
#' @param gt A gtable
#'
#' @return TRUE if none of its widths and heights are relative ("null") units
#'
#' @keywords internal
has_fixed_size <- function(gt) {
    !any(grid::unitType(gt$widths) == "null") && !any(grid::unitType(gt$heights) == "null")
}


#' @title Save a plot to a file sized to fit it
#'
#' @description Saves a ggplot or patchwork to PDF (vector), PNG or any format
#'              \code{\link[ggplot2]{ggsave}} supports.
#'
#'              For plots whose panels have a fixed size (all plots of this
#'              package, or any ggplot with \code{\link{panel_size_theme}}),
#'              the width and height of the file are computed from the panels,
#'              strips, legends and margins, so that nothing is cropped and no
#'              white space is added. They are measured on the same kind of
#'              device the file is written with, so that font metrics match.
#'              Patchworks of such plots, including nested ones, keep the
#'              size of every panel.
#'
#'              Plots without fixed panel sizes are saved at `width` by
#'              `height` inches, as with \code{ggsave()}.
#'
#' @param plot A ggplot or patchwork object
#'
#' @param filename File to save to; the format is inferred from the extension
#'                 (".pdf", ".png", ".svg", ".tiff", ...)
#'
#' @param width,height Size of the file in inches, used only for plots without
#'                     fixed panel sizes, for which they are required. Default
#'                     is NULL
#'
#' @param dpi Resolution of raster formats, default is 1000
#'
#' @param bg Background colour. If NULL (default), taken from the plot. A
#'           transparent background is kept in PDF, PNG, TIFF and SVG; formats
#'           without transparency (JPEG, BMP) get a white one instead
#'
#' @param ... Additional arguments passed to \code{\link[ggplot2]{ggsave}}
#'
#' @return The width and height of the saved file in inches, invisibly
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point() +
#'     theme_bioToolkit() + panel_size_theme(2.5, 2)
#' f <- tempfile(fileext = ".pdf")
#' save_plot(p, f)
#'
#' @seealso \code{\link{panel_size_theme}}, \code{\link{export_source_data}}
#'
#' @export
#'
save_plot <- function(plot, filename, width = NULL, height = NULL, dpi = 1000, bg = NULL, ...) {
    if (!(is_ggplot(plot) || inherits(plot, "patchwork"))) {
        stop("`plot` must be a ggplot or patchwork object", call. = FALSE)
    }
    # Build the gtable once, on the same kind of device it will be saved to, and
    # draw that exact gtable so measured and saved layouts cannot differ
    raster_ext <- grepl("\\.(png|jpe?g|tiff?|bmp)$", filename, ignore.case = TRUE)
    old_dev <- grDevices::dev.cur()
    if (raster_ext && requireNamespace("ragg", quietly = TRUE)) {
        tmp <- tempfile(fileext = ".png")
        ragg::agg_png(tmp, res = dpi)
        on.exit(unlink(tmp), add = TRUE)
    } else {
        grDevices::pdf(NULL)
    }
    measure_dev <- grDevices::dev.cur()
    is_patchwork <- inherits(plot, "patchwork")
    if (is_patchwork) plot <- patchwork_background(plot)
    # The gtable is measured while its device is still open; the device is
    # closed whatever happens, and the caller's device made current again
    measured <- tryCatch({
        gt <- if (is_patchwork) patchwork::patchworkGrob(fix_nested_patchworks(plot)) else ggplotGrob(plot)
        fixed <- has_fixed_size(gt)
        list(gt = gt, fixed = fixed,
             dims = if (fixed) fixed_gtable_size(gt, patchwork = is_patchwork))
    }, finally = {
        grDevices::dev.off(measure_dev)
        if (old_dev > 1) grDevices::dev.set(old_dev)
    })
    gt <- measured$gt
    if (measured$fixed) {
        dims <- measured$dims
        if (!is.null(width) || !is.null(height)) {
            message("The plot has fixed panel sizes; `width` and `height` are ignored")
        }
    } else {
        if (is.null(width) || is.null(height)) {
            stop("The plot has no fixed panel size: give `width` and `height` in inches, ",
                 "or add panel_size_theme() to the plot", call. = FALSE)
        }
        dims <- c(width = width, height = height)
    }

    if (is.null(bg)) {
        background <- gt$grobs[[match("background", gt$layout$name)]]
        bg <- if (is.null(background$gp$fill)) "transparent" else background$gp$fill
    }
    # JPEG and BMP have no alpha channel: a transparent background would be
    # written as black by some devices, so white is used instead
    if (grepl("\\.(jpe?g|bmp)$", filename, ignore.case = TRUE) &&
        (is.na(bg) || identical(bg, "transparent"))) {
        bg <- "white"
    }

    ggsave(filename, gt, width = dims[["width"]], height = dims[["height"]], units = "in",
           dpi = dpi, bg = bg, limitsize = FALSE, ...)
    invisible(dims)
}


#' Give a patchwork the background of its plots
#'
#' A patchwork draws its own background, taken from ggplot2's global theme
#' (white by default) unless set with \code{patchwork::plot_annotation()}. When
#' it is not set, the background of the first plot is used, so that a
#' patchwork of transparent plots stays transparent and one of white plots
#' stays white.
#'
#' @param plot A patchwork object
#'
#' @return The patchwork, with a background set
#'
#' @keywords internal
patchwork_background <- function(plot) {
    if (!is.null(plot$patches$annotation$theme$plot.background)) return(plot)
    first <- patchwork_parts(plot)$plots[[1]]
    if (inherits(first, "patchwork")) first <- patchwork_parts(first)$plots[[1]]
    if (!is_ggplot(first)) return(plot)
    background <- calc_element("plot.background", theme_get() + first$theme)
    if (is.null(background)) return(plot)
    plot + patchwork::plot_annotation(theme = theme(plot.background = background))
}


#' Give the nested patchworks of a patchwork a fixed size
#'
#' patchwork aligns the panels of the plots it combines, but a nested patchwork
#' (e.g. `a | b` in `(a | b) / c`) has no panel of its own, so its fixed size is
#' lost: the outer layout is left with relative sizes, and plots sharing a row
#' or column with it are stretched. On every level that holds a nested
#' patchwork, each plot is therefore drawn into a ggplot whose single panel has
#' exactly that plot's size, centred in its cell. Levels without nesting are
#' left to patchwork's own alignment. Must be called with the device the plot
#' is measured on open.
#'
#' @param plot A patchwork object
#'
#' @return The patchwork, rebuilt if it holds nested patchworks
#'
#' @keywords internal
fix_nested_patchworks <- function(plot) {
    parts <- patchwork_parts(plot)
    if (!any(vapply(parts$plots, inherits, logical(1), "patchwork"))) return(plot)
    plots <- lapply(parts$plots, function(p) {
        # Elements without a fixed size, such as spacers, are left as they are
        tryCatch(fixed_size_figure(p), error = function(e) p)
    })
    patchwork::wrap_plots(plots) + parts$layout + parts$annotation
}



#' The plots, layout and annotation of a patchwork
#'
#' A patchwork is the last plot added to it, carrying the earlier plots in its
#' `patches`; the last plot is split off and appended to them, unless it is the
#' empty placeholder that \code{patchwork::wrap_plots()} creates.
#'
#' @param plot A patchwork object
#'
#' @return A list with `plots`, `layout` and `annotation`
#'
#' @keywords internal
patchwork_parts <- function(plot) {
    parts <- plot$patches
    last <- plot
    last$patches <- NULL
    class(last) <- setdiff(class(last), "patchwork")
    if (!inherits(last, "plot_filler")) parts$plots <- c(parts$plots, list(last))
    parts
}


#' A ggplot holding a plot as one panel of the plot's fixed size
#'
#' @param plot A ggplot or patchwork object with fixed panel sizes
#'
#' @return A ggplot object
#'
#' @keywords internal
fixed_size_figure <- function(plot) {
    is_patchwork <- inherits(plot, "patchwork")
    gt <- if (is_patchwork) patchwork::patchworkGrob(fix_nested_patchworks(plot)) else ggplotGrob(plot)
    size <- fixed_gtable_size(gt, patchwork = is_patchwork)
    ggplot(data.frame(x = c(0, 1), y = c(0, 1)), aes(x = x, y = y)) +
        geom_blank() +
        annotation_custom(gt, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
        scale_x_continuous(expand = expansion(0)) +
        scale_y_continuous(expand = expansion(0)) +
        theme_void() +
        panel_size_theme(size[["width"]], size[["height"]]) +
        theme(plot.margin = margin(0, 0, 0, 0))
}


#' Absolute size of a gtable, including legends that overflow their cell
#'
#' A legend taller than a fixed-height panel (right/left) or wider than a
#' fixed-width panel (top/bottom) spills outside its cell and is not reflected
#' in the gtable's widths/heights. The legend is centred on its cell and the
#' gtable is centred on the device, so the device is enlarged to the legend size
#' plus twice the larger of the spaces before and after its cell, which keeps at
#' least the plot margin on both sides of the legend.
#'
#' In a patchwork, the space before and after a legend's cell holds the other
#' plots rather than margins, so the device is instead enlarged just enough for
#' the legend to keep `margin` inches from the device edge.
#'
#' @param gt A gtable with absolute widths and heights
#'
#' @param patchwork Whether `gt` is the gtable of a patchwork
#'
#' @param margin Smallest gap between a legend and the device edge in a
#'               patchwork, in inches
#'
#' @return Named numeric vector with width and height in inches
#'
#' @keywords internal
fixed_gtable_size <- function(gt, patchwork = FALSE, margin = 0.1) {
    if (any(grid::unitType(gt$widths) == "null") || any(grid::unitType(gt$heights) == "null")) {
        stop("The plot has no fixed size; set theme(panel.widths, panel.heights)", call. = FALSE)
    }
    widths <- grid::convertWidth(gt$widths, "in", valueOnly = TRUE)
    heights <- grid::convertHeight(gt$heights, "in", valueOnly = TRUE)
    width <- sum(widths)
    height <- sum(heights)

    for (i in grep("^guide-box", gt$layout$name)) {
        legend <- gt$grobs[[i]]
        if (inherits(legend, "zeroGrob")) next
        pos <- gt$layout[i, ]
        legend_w <- grid::convertWidth(grid::grobWidth(legend), "in", valueOnly = TRUE)
        legend_h <- grid::convertHeight(grid::grobHeight(legend), "in", valueOnly = TRUE)
        before_w <- sum(widths[seq_len(pos$l - 1)])
        after_w <- sum(widths[-seq_len(pos$r)])
        before_h <- sum(heights[seq_len(pos$t - 1)])
        after_h <- sum(heights[-seq_len(pos$b)])
        cell_w <- sum(widths[pos$l:pos$r])
        cell_h <- sum(heights[pos$t:pos$b])
        if (legend_w > cell_w + 1e-3) {
            width <- max(width, if (patchwork) {
                overflow_fit(legend_w, before_w + cell_w / 2, sum(widths), margin)
            } else {
                legend_w + 2 * max(before_w, after_w)
            })
        }
        if (legend_h > cell_h + 1e-3) {
            height <- max(height, if (patchwork) {
                overflow_fit(legend_h, before_h + cell_h / 2, sum(heights), margin)
            } else {
                legend_h + 2 * max(before_h, after_h)
            })
        }
    }
    c(width = width, height = height)
}


#' Device size fitting a legend centred on a point of a centred gtable
#'
#' @param legend Size of the legend
#'
#' @param centre Position of the legend's centre within the gtable
#'
#' @param total Size of the gtable
#'
#' @param margin Smallest gap between the legend and the device edge
#'
#' @return The device size
#'
#' @keywords internal
overflow_fit <- function(legend, centre, total, margin) {
    spill <- max(legend / 2 + margin - centre, centre + legend / 2 + margin - total, 0)
    total + 2 * spill
}
