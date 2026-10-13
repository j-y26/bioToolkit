# Internal dot-grid engine (roadmap item G0)
#
# Design
# ------
# The marker dot plot (`DotPlot_bioToolkit()`, G1), the DE dot plot
# (`plot_de_dots()`, G2) and the enrichment dot plot
# (`plot_enrichment_grid()`, G3) are the same plot: a grid of rows x columns
# of dots, each with a size and a fill, rows and columns optionally split into
# groups. Each plot function computes its own values and orders the rows and
# columns; `dot_grid_plot()` only draws them, so that the three look and size
# alike.
#
# Geometry (F2). Every row and column is a fixed number of inches
# (`cell.width`, `cell.height`). Discrete scales are expanded by exactly half
# a category (`grid_expansion()`), so that n categories fill n cells. Numbers
# outside 1..n on a discrete scale would widen its range and break this, so
# everything drawn outside the cells uses either its own facet or the panel
# edge (`-Inf` / `Inf`, which do not train scales):
#   * Groups are facets (`facet_grid(space = "free")`), one panel per group,
#     with the group name in the strip (top for column groups, left for row
#     groups). Column groups get a thin bracket line along the top panel
#     edge, drawn at +Inf in a gap added to the scale expansion
#     (`bracket.gap`); row groups are set apart by the breaks of the y-axis
#     line between their panels.
#   * The group markers (e.g. cell-type colours next to the cell-type rows)
#     are a facet of its own, with no axis label or tick, between the labels
#     and the axis: either one narrow column (or row) of tiles (a colour
#     strip), or one cell-wide column (or row) of coloured dots.
# Axes (2026-10-13). 1 pt black axis lines and ticks, as the "axis" plots,
# along the left and bottom edges of the data panels. ggplot2 draws axes
# around the outermost panels, which would put them outside the marker facet
# (labels, ticks, line, then the markers); so the theme axis lines and ticks
# are blank, the theme axis only places the labels, and the line and ticks
# are drawn by a layer (`GeomGridAxis`) on the data panels, with the markers
# outside them. Grid lines are likewise a layer drawn on the data panels
# only, so that none runs through the markers. The axes are set off from the
# first row and column of dots by a gap added to the expansion (`axis.gap`).
# Significance (`outline`): significant dots are drawn opaque with a dark
# ring, the others faded with a light ring, so that the two read apart at
# any dot size and the fill (e.g. fold change) of every dot is still shown.
# Both are explained by one legend (ring colour and opacity share its title
# and labels, so ggplot2 merges them).
# Dot sizes are areas (`scale_size_area()`), with the largest dot filling 85%
# of the smaller cell side, so that dots never touch and the size of a dot
# for a given value is the same in every figure with the same cell size.


#' Draw a grid of dots
#'
#' @param df One row per dot, with columns `x` and `y` (factors; their levels
#'           give the order of the columns from left to right and of the rows
#'           from top to bottom, and must be unique across groups), `size`
#'           (NA: no dot, but the row and column are kept),
#'           `fill`, and optionally `x_group` / `y_group` (factors: column /
#'           row groups), `outline` (logical: significant dots, drawn opaque
#'           with a dark ring; the others are faded) and `strip` (the colour
#'           of the row or column, for the group markers)
#'
#' @param x.labels,y.labels Display labels named by level of `x` / `y`
#'
#' @param x.italic,y.italic Whether each label is set in italics (genes):
#'                          logical vectors named by level, or a single value
#'
#' @param size.limits Range of `size` mapped to dot areas; values beyond are
#'                    drawn at the nearest end
#'
#' @param size.breaks,size.labels Breaks and labels of the size legend
#'
#' @param size.title Title of the size legend
#'
#' @param fill.colours Colours of the fill scale, from low to high
#'
#' @param fill.limits Limits of the fill scale; values beyond get the colour of
#'                    the nearest end
#'
#' @param fill.center Centre of a diverging fill scale, or NULL for a
#'                    sequential one
#'
#' @param fill.title Title of the fill legend
#'
#' @param na.value Fill of dots whose `fill` is NA
#'
#' @param outline.labels Legend labels of significant and other dots, when
#'                       `df$outline` is given
#'
#' @param outline.title Title of the significance legend
#'
#' @param strip.axis "x" when the group markers run along the columns (one
#'                   per column, below the grid), "y" when they run along the
#'                   rows (one per row, left of the grid), or NULL for none
#'
#' @param strip.style "strip" (a narrow band of tiles) or "dot" (one coloured
#'                    dot per row or column)
#'
#' @param cell.width,cell.height Size of one column / row in inches
#'
#' @param legend.position Position of the legends
#'
#' @param bg Background colour
#'
#' @return A ggplot object
#'
#' @keywords internal
dot_grid_plot <- function(df,
                          x.labels,
                          y.labels,
                          x.italic = FALSE,
                          y.italic = FALSE,
                          size.limits,
                          size.breaks,
                          size.labels = waiver(),
                          size.title,
                          fill.colours,
                          fill.limits,
                          fill.center = NULL,
                          fill.title,
                          na.value = "grey85",
                          outline.labels = NULL,
                          outline.title = NULL,
                          strip.axis = NULL,
                          strip.style = "strip",
                          cell.width = 0.18,
                          cell.height = 0.18,
                          legend.position = "right",
                          bg = "transparent") {
    check_fixed_panel_support()
    sz <- theme_sizes()
    lw <- pt_to_linewidth(sz$line)
    has_x_group <- !is.null(df$x_group)
    has_y_group <- !is.null(df$y_group)
    has_outline <- !is.null(df$outline)
    if (!is.null(strip.axis)) strip.axis <- match.arg(strip.axis, c("x", "y"))
    strip.style <- match.arg(strip.style, c("strip", "dot"))

    # Rows are listed top to bottom; ggplot2 draws the first level at the
    # bottom, so they are reversed
    df$x <- droplevels(df$x)
    df$y <- factor(df$y, levels = rev(levels(droplevels(df$y))))
    if (!has_outline) df$outline <- NULL
    if (!has_x_group) df$x_group <- factor("")
    if (!has_y_group) df$y_group <- factor("")
    df$x_group <- droplevels(df$x_group)
    df$y_group <- droplevels(df$y_group)

    # Gaps in inches: between the panels of groups; between the top panel
    # edge and the column-group bracket line; and between the axes (left and
    # bottom edges) and the first column / row of dots, so that the largest
    # dots do not touch the axis lines. The bracket and axis gaps are added to
    # the scale expansion of that edge, so they are part of the panel and the
    # cells keep their size
    spacing <- 0.08
    bracket.gap <- 0.06
    axis.gap <- 0.04
    x_bracket <- has_x_group
    # Axis line and ticks (pt), as the "axis" theme: 1 pt, ticks a third of
    # the tick-label size long
    axis.width <- sz$axis
    tick.length <- sz$small / 3

    # The colour strip is a facet of its own: an extra first column group
    # (strip along the rows) or last row group (strip along the columns)
    strip_level <- ".strip"
    strip_df <- NULL
    if (!is.null(strip.axis)) {
        if (is.null(df$strip)) stop("`strip.axis` needs a `strip` column", call. = FALSE)
        if (strip.axis == "y") {
            strip_df <- unique(df[c("y", "y_group", "strip")])
            strip_df$x <- strip_level
            strip_df$x_group <- strip_level
        } else {
            strip_df <- unique(df[c("x", "x_group", "strip")])
            strip_df$y <- strip_level
            strip_df$y_group <- strip_level
        }
    }
    x_lv <- levels(df$x)
    y_lv <- levels(df$y)
    xg_lv <- levels(df$x_group)
    yg_lv <- levels(df$y_group)
    if (!is.null(strip.axis)) {
        if (strip.axis == "y") {
            x_lv <- c(strip_level, x_lv)
            xg_lv <- c(strip_level, xg_lv)
        } else {
            # Bottom row: first level, since rows are drawn bottom-up; its
            # group is the last facet row
            y_lv <- c(strip_level, y_lv)
            yg_lv <- c(yg_lv, strip_level)
        }
    }
    relevel_all <- function(d) {
        d$x <- factor(as.character(d$x), levels = x_lv)
        d$y <- factor(as.character(d$y), levels = y_lv)
        d$x_group <- factor(as.character(d$x_group), levels = xg_lv)
        d$y_group <- factor(as.character(d$y_group), levels = yg_lv)
        d
    }
    df <- relevel_all(df)
    if (!is.null(strip_df)) strip_df <- relevel_all(strip_df)

    # Panel sizes: number of categories per group times the cell size; the
    # strip facet is a narrow band (or one cell, for dots); panels next to a
    # bracket get its gap
    strip_size <- if (strip.style == "dot") {
        if (identical(strip.axis, "y")) cell.width else cell.height
    } else {
        0.1
    }
    n_x <- vapply(xg_lv, function(g) length(unique(df$x[df$x_group == g])), integer(1))
    n_y <- vapply(yg_lv, function(g) length(unique(df$y[df$y_group == g])), integer(1))
    # Gaps in categories (the unit of the expansion)
    gap_x <- axis.gap / cell.width
    gap_bottom <- axis.gap / cell.height
    gap_top <- if (x_bracket) bracket.gap / cell.height else 0
    widths <- grid_panel_size(n_x, cell.width) + gap_x * cell.width
    heights <- grid_panel_size(n_y, cell.height) + (gap_bottom + gap_top) * cell.height
    # The strip facet holds one category plus the gaps, in a narrow band, so
    # that its tiles keep the strip's thickness (the expansion is shared by
    # all facets, so the strip facet has the gaps too, scaled to its size)
    if (!is.null(strip.axis)) {
        if (strip.axis == "y") {
            widths[1] <- strip_size * (1 + gap_x)
        } else {
            heights[length(heights)] <- strip_size * (1 + gap_bottom + gap_top)
        }
    }
    # Expansion: half a category on each side, plus the axis gap on the left
    # and bottom and the bracket gap on the top
    x_expand <- expansion(add = c(0.5 + gap_x, 0.5))
    y_expand <- expansion(add = c(0.5 + gap_bottom, 0.5 + gap_top))

    # Largest dot: 85% of the smaller cell side. A filled circle (shape 21)
    # of size s is 0.75 * s * .pt points across
    max_size <- 0.85 * min(cell.width, cell.height) * 72 / (0.75 * ggplot2::.pt)

    label_fun <- function(labels, italic) {
        function(b) {
            b <- b[b != strip_level]
            lab <- unname(labels[b])
            it <- if (length(italic) == 1 && is.null(names(italic))) rep(italic, length(b)) else italic[b]
            italic_labels(lab, unname(!is.na(it) & it))
        }
    }
    no_strip <- function(lim) lim[lim != strip_level]

    # Column-group labels that are wider than their panel are turned upright,
    # so that neighbouring labels never overlap. Widths are estimated from the
    # number of characters (bold text, about 0.6 em per character)
    x_group_angle <- 0
    if (has_x_group) {
        real <- xg_lv != strip_level
        label_w <- nchar(xg_lv[real]) * 0.6 * sz$small / 72
        if (any(label_w > widths[real] + spacing)) x_group_angle <- getOption("bt.angle", 90)
    }

    # Every row and column is kept, even where no dot is drawn (`size` NA,
    # e.g. below `pct.min`): free facet scales would otherwise drop it
    grid_df <- df
    dots <- df[!is.na(df$size), , drop = FALSE]
    p <- ggplot(dots, aes(x = .data$x, y = .data$y)) +
        geom_blank(data = grid_df)
    # Grid lines, one per row and column, on the data panels only: each line
    # is given the facets (column and row group) of the panels it crosses, so
    # that none is drawn through the marker facet
    v_lines <- unique(grid_df[c("x", "x_group", "y_group")])
    h_lines <- unique(grid_df[c("y", "x_group", "y_group")])
    p <- p +
        geom_segment(aes(x = .data$x, xend = .data$x, y = -Inf, yend = Inf), data = v_lines,
                     inherit.aes = FALSE, colour = "grey92", linewidth = lw) +
        geom_segment(aes(x = -Inf, xend = Inf, y = .data$y, yend = .data$y), data = h_lines,
                     inherit.aes = FALSE, colour = "grey92", linewidth = lw)
    if (!is.null(strip_df)) {
        if (strip.style == "strip") {
            p <- p + geom_tile(aes(fill = NULL), data = strip_df, fill = strip_df$strip, colour = NA,
                               width = 1, height = 1)
        } else {
            # A dot a little smaller than the largest data dot, so that it reads
            # as a label of the row (column) rather than as a value
            # (shape 16 with a set colour: a point with an NA colour is dropped)
            p <- p + geom_point(aes(fill = NULL, size = NULL), data = strip_df, colour = strip_df$strip,
                                shape = 16, size = max_size * 0.75)
        }
    }
    if (has_outline) {
        # Faded dots first, then significant ones on top. The ring colour and
        # opacity share one legend (same title and labels)
        ring <- c(pt_to_stroke(0.7), pt_to_stroke(0.3))
        ring_cols <- stats::setNames(c("black", "grey45"), outline.labels)
        ring_alpha <- stats::setNames(c(1, 0.3), outline.labels)
        plain <- dots[!dots$outline, , drop = FALSE]
        ringed <- dots[dots$outline, , drop = FALSE]
        plain$ring <- factor(rep(outline.labels[2], nrow(plain)), levels = outline.labels)
        ringed$ring <- factor(rep(outline.labels[1], nrow(ringed)), levels = outline.labels)
        sig_guide <- guide_legend(order = 3, override.aes = list(
            size = max_size * 0.7, fill = "grey50", stroke = ring
        ))
        p <- p +
            geom_point(aes(size = .data$size, fill = .data$fill, colour = .data$ring, alpha = .data$ring),
                       data = plain, shape = 21, stroke = ring[2]) +
            geom_point(aes(size = .data$size, fill = .data$fill, colour = .data$ring, alpha = .data$ring),
                       data = ringed, shape = 21, stroke = ring[1]) +
            scale_colour_manual(values = ring_cols, drop = FALSE, name = outline.title, guide = sig_guide) +
            # Merged with the colour legend, whose key overrides it carries
            scale_alpha_manual(values = ring_alpha, drop = FALSE, name = outline.title,
                               guide = guide_legend(order = 3))
    } else {
        p <- p + geom_point(aes(size = .data$size, fill = .data$fill), shape = 21,
                            colour = "grey40", stroke = pt_to_stroke(0.25))
    }

    # Gap (pt) between the tick labels and the ticks, or the markers when
    # they sit between the labels and the axis
    label_gap <- function(markers) if (markers) sz$small / 4 else tick.length + sz$small / 4

    vertical <- !(legend.position %in% c("top", "bottom"))
    p <- p +
        scale_x_discrete(labels = label_fun(x.labels, x.italic), breaks = no_strip,
                         expand = x_expand, drop = TRUE) +
        scale_y_discrete(labels = label_fun(y.labels, y.italic), breaks = no_strip,
                         expand = y_expand, drop = TRUE) +
        scale_size_area(
            max_size = max_size, limits = size.limits, breaks = size.breaks,
            labels = size.labels, oob = scales::oob_squish, name = size.title,
            guide = guide_legend(order = 2, override.aes = list(fill = "grey50", colour = "grey25"))
        ) +
        scale_fill_gradientn(
            colours = fill.colours, limits = fill.limits, oob = scales::oob_squish,
            rescaler = if (is.null(fill.center)) {
                function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) scales::rescale(x, to, from)
            } else {
                function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
                    scales::rescale_mid(x, to, from, mid = fill.center)
                }
            },
            na.value = na.value, name = fill.title,
            guide = guide_colourbar(
                order = 1,
                direction = if (vertical) "vertical" else "horizontal",
                theme = theme(
                    legend.key.height = unit(if (vertical) 0.9 else 0.12, "in"),
                    legend.key.width = unit(if (vertical) 0.12 else 0.9, "in"),
                    legend.frame = element_rect(colour = "black", linewidth = 0.3),
                    legend.ticks = element_line(colour = "black", linewidth = 0.3),
                    legend.ticks.length = unit(0.03, "in")
                )
            )
        ) +
        theme_bioToolkit("grid", bg = bg) +
        panel_size_theme(widths, heights, panel.spacing = spacing) +
        # The axis lines and ticks hang outside the panels (layer below)
        coord_cartesian(clip = "off") +
        theme(
            # Axis lines, ticks and grid lines are layers on the data panels
            # (see the notes at the top); the theme axes only place the
            # labels: past the ticks, or next to the markers on their side
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.ticks.length = unit(0, "pt"),
            panel.grid.major = element_blank(),
            axis.text.x = element_text(margin = margin(t = label_gap(identical(strip.axis, "x")))),
            axis.text.y = element_text(margin = margin(r = label_gap(identical(strip.axis, "y")))),
            legend.position = legend.position,
            legend.title = element_text(size = sz$small, margin = margin(b = 0.15, unit = "cm")),
            legend.title.position = "top",
            legend.key.size = unit(max(min(cell.width, cell.height), 0.16), "in"),
            legend.spacing.y = unit(0.15, "in"),
            # Groups are named by their strips, set off from the grid by the
            # bracket lines
            strip.placement = "outside",
            strip.clip = "off",
            strip.text.x.top = element_text(size = sz$small, face = "bold", angle = x_group_angle,
                                            hjust = if (x_group_angle == 0) 0.5 else 0,
                                            vjust = if (x_group_angle == 0) 0.5 else 0.5,
                                            margin = margin(b = 2, t = 1)),
            strip.text.y.left = element_text(size = sz$small, face = "bold", angle = 0, hjust = 1,
                                             margin = margin(r = 3, l = 1))
        )

    if (has_x_group || has_y_group || !is.null(strip.axis)) {
        p <- p + facet_grid(
            rows = vars(.data$y_group), cols = vars(.data$x_group),
            scales = "free", space = "free",
            switch = "y",
            labeller = as_labeller(function(x) ifelse(x == strip_level, "", x))
        )
        # Facets that exist only to hold the strip have no label
        if (!has_x_group) p <- p + theme(strip.text.x = element_blank())
        if (!has_y_group) p <- p + theme(strip.text.y = element_blank())
    }

    # Bracket lines above column groups (top row of panels), drawn at +Inf so
    # they do not train the scales
    data_xg <- setdiff(xg_lv, strip_level)
    data_yg <- setdiff(yg_lv, strip_level)
    if (has_x_group) {
        br <- data.frame(x_group = factor(data_xg, levels = xg_lv),
                         y_group = factor(yg_lv[1], levels = yg_lv),
                         x = -Inf, xend = Inf, y = Inf, yend = Inf)
        p <- p + geom_segment(data = br, aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
                              inherit.aes = FALSE, colour = "black", linewidth = pt_to_linewidth(0.75))
    }

    # Axes: the y-axis on the first column of data panels (one tick per row),
    # the x-axis on the last row of data panels (one tick per column). The
    # y-axis line stops below the bracket gap, so that it does not run into
    # the column-group brackets
    y_ticks <- unique(grid_df[c("y", "y_group")])
    y_ticks$x_group <- factor(data_xg[1], levels = xg_lv)
    y_ticks$x <- -Inf
    x_ticks <- unique(grid_df[c("x", "x_group")])
    x_ticks$y_group <- factor(data_yg[length(data_yg)], levels = yg_lv)
    x_ticks$y <- -Inf
    axis_params <- list(width = axis.width, length = tick.length)
    p + grid_axis_layer(y_ticks, c(axis_params, side = "left", top.trim = gap_top * cell.height)) +
        grid_axis_layer(x_ticks, c(axis_params, side = "bottom", top.trim = 0))
}


# Axis line and ticks of the dot grid, drawn by a layer on the data panels
# (see the notes at the top). One row of data per tick: `y` for the left axis
# or `x` for the bottom axis (the other position is -Inf), and the facets of
# the panel. Panels without a row get no axis
GeomGridAxis <- ggproto("GeomGridAxis", Geom,
    required_aes = c("x", "y"),
    draw_panel = function(data, panel_params, coord, side = "left", width = 1, length = 3, top.trim = 0) {
        coords <- coord$transform(data, panel_params)
        # Widths in points; grid line widths are in 1/96 inch
        lwd <- width * 96 / 72
        zero <- unit(0, "npc")
        out <- unit(0, "npc") - unit(length, "pt")
        if (side == "left") {
            at <- unit(coords$y, "npc")
            ticks <- grid::segmentsGrob(zero, at, out, at)
            line <- grid::segmentsGrob(zero, zero, zero, unit(1, "npc") - unit(top.trim, "in"))
        } else {
            at <- unit(coords$x, "npc")
            ticks <- grid::segmentsGrob(at, zero, at, out)
            line <- grid::segmentsGrob(zero, zero, unit(1, "npc"), zero)
        }
        # Butt ends on the ticks; square ends on the line, so that the two
        # axes meet in a full corner
        ticks$gp <- grid::gpar(col = "black", lwd = lwd, lineend = "butt")
        line$gp <- grid::gpar(col = "black", lwd = lwd, lineend = "square")
        grid::grobTree(ticks, line)
    }
)


#' Layer of the dot-grid axes
#'
#' @param data One row per tick (see `GeomGridAxis`)
#'
#' @param params `side` ("left" or "bottom"), `width` and `length` of the
#'               line and ticks (pt), and `top.trim` (inches left out at the
#'               top of a left axis line)
#'
#' @return A ggplot2 layer
#'
#' @keywords internal
grid_axis_layer <- function(data,
                            params) {
    layer(geom = GeomGridAxis, stat = "identity", position = "identity", data = data,
          mapping = aes(x = .data$x, y = .data$y), inherit.aes = FALSE, show.legend = FALSE,
          params = params)
}


#' Order the rows of a matrix by hierarchical clustering, within groups
#'
#' Rows are clustered with Ward's method (`ward.D2`) on Euclidean distances,
#' separately within each group, and groups keep their order. Missing values
#' count as 0. Dendrograms are not drawn; only the order is used.
#'
#' @param mat A numeric matrix with row names
#'
#' @param groups A factor along the rows, or NULL for one group
#'
#' @return The row names, in clustered order
#'
#' @keywords internal
cluster_order <- function(mat,
                          groups = NULL) {
    mat[!is.finite(mat)] <- 0
    if (is.null(groups)) groups <- factor(rep("", nrow(mat)))
    groups <- droplevels(as.factor(groups))
    unlist(lapply(levels(groups), function(g) {
        rows <- which(groups == g)
        if (length(rows) < 3) return(rownames(mat)[rows])
        hc <- stats::hclust(stats::dist(mat[rows, , drop = FALSE]), method = "ward.D2")
        rownames(mat)[rows][hc$order]
    }), use.names = FALSE)
}
