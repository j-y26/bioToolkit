# Legends as separate figures
#
# Design
# ------
# * Manual figure layout (e.g. in Illustrator) often needs the legend on its
#   own: `extract_legend()` takes it out of a plot, and `save_plot()` saves it
#   to a file sized to fit, like a plot.
# * By default the legend is exactly as it is drawn in the plot. Changes
#   (direction, title position, text sizes, ...) are made by re-drawing the
#   legend with ggplot2's own theme settings before taking it out, never by
#   editing the drawn legend, so anything ggplot2 can do to a legend works.
# * ggplot2 measures text when it builds a plot, with the font metrics of the
#   open device; a legend built for a PDF and written to a PNG (different
#   fonts) gets clipped text. The extracted legend therefore keeps the plot it
#   came from, and `save_plot()` re-builds it on the kind of device it saves
#   with, as it does for plots.
# * The legend is returned as a grob, so it can also be drawn with
#   `grid::grid.draw()` or placed in a patchwork.


#' @title Extract the legend of a plot
#'
#' @description Takes the legend out of a ggplot, to save it as its own file
#'              with \code{\link{save_plot}} or place it in a layout. By default
#'              the legend is exactly as drawn in the plot; the arguments
#'              change it when needed.
#'
#' @param plot A ggplot
#'
#' @param direction Direction of the legend: NULL (default) keeps the plot's;
#'                  "horizontal" lays the keys out in a row, "vertical" in a
#'                  column. Several legends (e.g. colour and size) are placed
#'                  side by side or stacked accordingly
#'
#' @param title.position Position of the legend title: NULL (default) keeps
#'                       the plot's; or one of "top", "bottom", "left" and
#'                       "right"
#'
#' @param ... Other theme settings applied to the legend before it is taken
#'            out, e.g. `legend.text = element_text(size = 8)`,
#'            `legend.key.spacing.x = unit(4, "pt")` or `legend.box = "vertical"`
#'
#' @details
#' A plot whose legend is hidden (`legend.position = "none"`, as in
#' \code{\link{plot_volcano}} with `annotate.counts = "arrow"`) still gives its
#' legend, drawn as it would be on the right. Legends of patchworks are not
#' supported; extract them from the individual plots.
#'
#' The legend saved by \code{\link{save_plot}} is re-built on the kind of device
#' it is written with, so its text fits whatever the file format. A small
#' margin (the theme's small text size, 9 pt) is kept around it, so that the end
#' labels of colour bars, which overhang the legend box, are not cut.
#'
#' In a patchwork the legend has no fixed size: give `width` and `height` to
#' \code{\link{save_plot}}, or save the plot (without its legend,
#' `plot + theme(legend.position = "none")`) and the legend separately.
#'
#' @return The legend as a grob (class `bt_legend`), drawn when printed
#'
#' @examples
#' set.seed(1)
#' de <- data.frame(
#'     logFC = rnorm(500, sd = 1.5),
#'     PValue = runif(500)^4,
#'     row.names = paste0("Gene", 1:500)
#' )
#' de$FDR <- p.adjust(de$PValue, "BH")
#' p <- plot_volcano(de, method = "edgeR", label = "none")
#' leg <- extract_legend(p)
#' leg
#' extract_legend(p, direction = "horizontal")
#' \dontrun{
#' save_plot(p + ggplot2::theme(legend.position = "none"), "volcano.pdf")
#' save_plot(leg, "volcano_legend.pdf")
#' }
#'
#' @seealso \code{\link{save_plot}}
#'
#' @export
#'
extract_legend <- function(plot,
                           direction = NULL,
                           title.position = NULL,
                           ...) {
    if (!is_ggplot(plot) || inherits(plot, "patchwork")) {
        stop("`plot` must be a ggplot (not a patchwork)", call. = FALSE)
    }
    if (!is.null(direction)) {
        direction <- match.arg(direction, c("vertical", "horizontal"))
        # legend.box arranges several legends the same way as the keys
        plot <- plot + theme(legend.direction = direction, legend.box = direction)
    }
    if (!is.null(title.position)) {
        title.position <- match.arg(title.position, c("top", "bottom", "left", "right"))
        plot <- plot + theme(legend.title.position = title.position)
    }
    if (...length()) plot <- plot + theme(...)

    # Built once now (font metrics of a PDF device) so that the legend can be
    # drawn and placed in layouts; save_plot() re-builds it for its device
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    built <- legend_gtable(plot)
    legend <- built$legend
    attr(legend, "legend_plot") <- built$plot
    class(legend) <- c("bt_legend", class(legend))
    legend
}


#' Build a plot and take out its legend
#'
#' Builds the plot on the current device. A hidden legend is drawn on the
#' right instead.
#'
#' @param plot A ggplot
#'
#' @return A list with the `legend` (a gtable), the full plot gtable (`gt`)
#'         and the `plot` it was built from (with the legend shown)
#'
#' @keywords internal
legend_gtable <- function(plot) {
    find <- function(gt) {
        boxes <- grep("^guide-box", gt$layout$name)
        boxes[!vapply(gt$grobs[boxes], inherits, logical(1), "zeroGrob")]
    }
    gt <- ggplotGrob(plot)
    boxes <- find(gt)
    if (!length(boxes)) {
        # Hidden legend (legend.position = "none"): draw it on the right
        plot <- plot + theme(legend.position = "right")
        gt <- ggplotGrob(plot)
        boxes <- find(gt)
    }
    if (!length(boxes)) stop("The plot has no legend", call. = FALSE)
    if (length(boxes) > 1) {
        message("The plot has legends in ", length(boxes), " positions; the first (",
                sub("guide-box-", "", gt$layout$name[boxes[1]]), ") is extracted")
    }
    list(legend = gt$grobs[[boxes[1]]], gt = gt, plot = plot)
}


#' Size of a legend in inches, with a margin
#'
#' Measured on the current device.
#'
#' @param legend A legend gtable
#'
#' @return `c(width, height)` in inches, including a margin of the theme's
#'         small text size (9 pt) on each side
#'
#' @keywords internal
legend_size <- function(legend) {
    # Colour bar end labels are centred on the bar ends and overhang the
    # legend box; the margin keeps them inside the file
    pad <- theme_sizes()$small / 72
    c(width = grid::convertWidth(grid::grobWidth(legend), "in", valueOnly = TRUE) + 2 * pad,
      height = grid::convertHeight(grid::grobHeight(legend), "in", valueOnly = TRUE) + 2 * pad)
}


#' @param x A legend from `extract_legend()`
#'
#' @rdname extract_legend
#' @export
print.bt_legend <- function(x,
                            ...) {
    grid::grid.newpage()
    grid::grid.draw(x)
    invisible(x)
}


# Legend blocks (revision of 2026-10-14)
# ---------------------------------------
# Legends of several plots arranged in columns, e.g. the colour bars of the
# panels of plot_enrichment_dots() in one column and their size legends in
# the next, each column in panel order. Like extract_legend(), the legends
# are built from their plots on the device they are drawn on (in
# makeContent()), so that their text fits whatever the file format; the
# size used for the layout is measured once on a PDF device, with some room
# to spare.


#' A block of legends in columns
#'
#' @param columns A list of columns, each a list of ggplots whose legends are
#'                stacked top to bottom
#'
#' @param gap Gap between legends and between columns, in inches
#'
#' @return A grob (class `bt_legend_block`), its legends built when drawn
#'
#' @keywords internal
legend_block_grob <- function(columns,
                              gap = 0.15) {
    columns <- Filter(length, columns)
    grid::gTree(columns = columns, gap = gap, cl = "bt_legend_block")
}


#' Legends of a legend block, built and placed on the current device
#'
#' @param x A legend block
#'
#' @return A list of `grobs` (legends in viewports, top-left aligned),
#'         `width` and `height` of the block, in inches
#'
#' @keywords internal
legend_block_layout <- function(x) {
    kids <- list()
    x0 <- 0
    height <- 0
    for (col in x$columns) {
        legends <- lapply(col, function(p) legend_gtable(p)$legend)
        ws <- vapply(legends, function(l) grid::convertWidth(grid::grobWidth(l), "in", valueOnly = TRUE),
                     numeric(1))
        hs <- vapply(legends, function(l) grid::convertHeight(grid::grobHeight(l), "in", valueOnly = TRUE),
                     numeric(1))
        y0 <- 0
        for (k in seq_along(legends)) {
            vp <- grid::viewport(
                x = grid::unit(x0, "in"),
                y = grid::unit(1, "npc") - grid::unit(y0, "in"),
                width = grid::unit(ws[k], "in"),
                height = grid::unit(hs[k], "in"),
                just = c("left", "top")
            )
            kids <- c(kids, list(grid::grobTree(legends[[k]], vp = vp)))
            y0 <- y0 + hs[k] + x$gap
        }
        height <- max(height, y0 - x$gap)
        x0 <- x0 + max(ws) + x$gap
    }
    list(grobs = kids, width = x0 - x$gap, height = height)
}


#' @importFrom grid makeContent
#' @export
makeContent.bt_legend_block <- function(x) {
    grid::setChildren(x, do.call(grid::gList, legend_block_layout(x)$grobs))
}


#' A ggplot of fixed width holding a legend block at its top left
#'
#' Only the width is fixed: placed in a patchwork over several rows of
#' fixed-size panels, the figure takes their height.
#'
#' @param block A legend block (\code{legend_block_grob()})
#'
#' @param pad Room left of the legends, in inches
#'
#' @return A ggplot object, with the block's measured height (inches, PDF
#'         metrics plus 5%) as attribute `block_height`
#'
#' @keywords internal
legend_block_figure <- function(block,
                                pad = 0.15) {
    # Measured on a PDF device; 5% more, so that other devices' fonts fit
    grDevices::pdf(NULL)
    size <- tryCatch(legend_block_layout(block), finally = grDevices::dev.off())
    placed <- grid::grobTree(block, vp = grid::viewport(
        x = grid::unit(pad, "in"),
        width = grid::unit(1, "npc") - grid::unit(pad, "in"),
        just = "left"
    ))
    fig <- ggplot(data.frame(x = c(0, 1), y = c(0, 1)), aes(x = .data$x, y = .data$y)) +
        geom_blank() +
        annotation_custom(placed, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
        scale_x_continuous(expand = expansion(0)) +
        scale_y_continuous(expand = expansion(0)) +
        coord_cartesian(clip = "off") +
        theme_void() +
        theme(
            panel.widths = grid::unit(size$width * 1.05 + pad, "in"),
            plot.margin = margin(0, 0, 0, 0)
        )
    attr(fig, "block_height") <- size$height * 1.05
    fig
}
