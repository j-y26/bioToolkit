# Group annotation bars (revision of 2026-10-14)
#
# Design
# ------
# A higher-level grouping of the rows or columns of a plot (genes grouped by
# cell type, `split.by` blocks, `groups` of comparisons) is drawn the same way
# in every plot of the package: a bar along the cells of the group, set a
# little off the grid, with the group's name in bold beyond it. Decisions
# (user, 2026-10-14):
#
# * The bar is 3.5 pt thick by default (`group.bar.width`): about a quarter of
#   a default 0.18 in (13 pt) cell, so it reads as a label band, not as an
#   axis line, and is wide enough to show a colour. It is black by default,
#   or coloured by group through the plot's `*.cols` argument (the same
#   values as `cols`: a qualitative palette name, or colours named by group).
# * By default the bar and name sit on the side opposite the axis labels
#   (e.g. genes labelled at the bottom: gene groups on top). `*.side =
#   "labels"` moves them to the side of the labels, beyond the labels.
#
# Geometry. The bar is drawn by a layer (`GeomGroupBar`) on the outermost
# panel of that side, outside the panel (the plot is drawn with
# `clip = "off"`): `GROUP_BAR_GAP` points from the panel edge, or from the far
# edge of the tick labels on the label side, whose size is measured when the
# plot is drawn (`label_extent()`). The group names are the facet strips;
# their margin on the panel side leaves room for the bar
# (`group_bar_margin()`). The bar spans exactly the cells of the group, so
# the bars of neighbouring groups are set apart by the gap between panels.


# Gaps (pt): between the grid (or the tick labels) and the bar, and between
# the bar and the group name
GROUP_BAR_GAP <- 2.5
GROUP_BAR_TEXT_GAP <- 2


# Bar along one side of a panel, outside it. One row of data per bar, with
# `x` / `xend` (bars on top or bottom) or `y` / `yend` (bars left or right):
# the positions of the first and last category of the group, and `fill`
GeomGroupBar <- ggproto("GeomGroupBar", Geom,
    default_aes = aes(fill = "black"),
    optional_aes = c("x", "xend", "y", "yend"),
    draw_key = draw_key_rect,
    draw_panel = function(data,
                          panel_params,
                          coord,
                          side = "top",
                          width = 3.5,
                          offset = grid::unit(GROUP_BAR_GAP, "pt"),
                          half = 0.5) {
        horizontal <- side %in% c("top", "bottom")
        # The bar spans its categories, each `half` a unit on either side of
        # its position
        if (horizontal) {
            lo <- pmin(data$x, data$xend) - half
            hi <- pmax(data$x, data$xend) + half
            ends <- coord$transform(data.frame(x = c(lo, hi), y = 0), panel_params)$x
        } else {
            lo <- pmin(data$y, data$yend) - half
            hi <- pmax(data$y, data$yend) + half
            ends <- coord$transform(data.frame(x = 0, y = c(lo, hi)), panel_params)$y
        }
        n <- nrow(data)
        from <- grid::unit(ends[seq_len(n)], "npc")
        to <- grid::unit(ends[n + seq_len(n)], "npc")
        thick <- grid::unit(width, "pt")
        gp <- grid::gpar(fill = data$fill, col = NA)
        switch(
            side,
            top = grid::rectGrob(x = from, width = to - from, y = grid::unit(1, "npc") + offset,
                                 height = thick, just = c("left", "bottom"), gp = gp),
            bottom = grid::rectGrob(x = from, width = to - from, y = grid::unit(0, "npc") - offset,
                                    height = thick, just = c("left", "top"), gp = gp),
            right = grid::rectGrob(y = from, height = to - from, x = grid::unit(1, "npc") + offset,
                                   width = thick, just = c("left", "bottom"), gp = gp),
            left = grid::rectGrob(y = from, height = to - from, x = grid::unit(0, "npc") - offset,
                                  width = thick, just = c("right", "bottom"), gp = gp)
        )
    }
)


#' Layer of group annotation bars
#'
#' @param data One row per bar: the facet columns of the panel it is drawn
#'             on, `x` / `xend` (bars on top or bottom) or `y` / `yend` (bars
#'             left or right) giving the first and last category of the group,
#'             and `fill` (colours)
#'
#' @param side Side of the panel: "top", "bottom", "left" or "right"
#'
#' @param width Thickness of the bar in points
#'
#' @param offset Distance from the panel edge to the bar, a grid unit
#'
#' @param half Half the width of a category, in data units (0.5 on a discrete
#'             scale)
#'
#' @return A ggplot2 layer
#'
#' @keywords internal
group_bar_layer <- function(data,
                            side,
                            width = 3.5,
                            offset = grid::unit(GROUP_BAR_GAP, "pt"),
                            half = 0.5) {
    side <- match.arg(side, c("top", "bottom", "left", "right"))
    mapping <- if (side %in% c("top", "bottom")) {
        aes(x = .data$x, xend = .data$xend, fill = I(.data$fill))
    } else {
        aes(y = .data$y, yend = .data$yend, fill = I(.data$fill))
    }
    layer(
        geom = GeomGroupBar,
        stat = "identity",
        position = "identity",
        data = data,
        mapping = mapping,
        inherit.aes = FALSE,
        show.legend = FALSE,
        params = list(side = side, width = width, offset = offset, half = half)
    )
}


#' Margin of a group name, leaving room for its bar
#'
#' @param side Side of the panel the name is on
#'
#' @param width Thickness of the bar in points
#'
#' @return A `margin()`, in points
#'
#' @keywords internal
group_bar_margin <- function(side,
                             width) {
    # Room for the bar between the panel (or labels) and the name
    room <- GROUP_BAR_GAP + width + GROUP_BAR_TEXT_GAP
    switch(
        side,
        top = margin(b = room, t = 1, unit = "pt"),
        bottom = margin(t = room, b = 1, unit = "pt"),
        left = margin(r = room, l = 1, unit = "pt"),
        right = margin(l = room, r = 1, unit = "pt")
    )
}


#' Size of a set of tick labels across the axis, measured when drawn
#'
#' @param labels Labels (strings or plotmath expressions)
#'
#' @param size Font size in points
#'
#' @param angle Rotation of the labels in degrees
#'
#' @param across "height" for labels of a horizontal axis (x), "width" for
#'               those of a vertical axis (y)
#'
#' @return A grid unit
#'
#' @keywords internal
label_extent <- function(labels,
                         size,
                         angle = 0,
                         across = c("height", "width")) {
    across <- match.arg(across)
    if (!length(labels)) return(grid::unit(0, "pt"))
    # The largest of the labels measured one by one (a text grob of several
    # rotated plotmath labels is not measured by its longest label)
    sizes <- lapply(seq_along(labels), function(i) {
        g <- grid::textGrob(labels[i], rot = angle, gp = grid::gpar(fontsize = size))
        if (across == "height") grid::unit(1, "grobheight", g) else grid::unit(1, "grobwidth", g)
    })
    max(do.call(grid::unit.c, sizes))
}


#' Colours of group bars
#'
#' @param cols NULL (default: black), a qualitative palette name, or colours
#'             named by group (or in level order)
#'
#' @param groups Group levels
#'
#' @param arg Name of the argument, for messages
#'
#' @return Colours named by group
#'
#' @keywords internal
group_bar_colours <- function(cols,
                              groups,
                              arg) {
    if (is.null(cols)) return(stats::setNames(rep("black", length(groups)), groups))
    group_colours(cols, factor(groups, levels = groups), arg = arg)$cols
}


#' Check the thickness of group bars
#'
#' @param width The thickness in points
#'
#' @return `width`, invisibly
#'
#' @keywords internal
check_bar_width <- function(width) {
    if (!(is.numeric(width) && length(width) == 1 && !is.na(width) && width > 0)) {
        stop("`group.bar.width` must be a single positive number (points)", call. = FALSE)
    }
    invisible(width)
}
