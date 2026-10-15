# Cell composition between two conditions (roadmap item S3)
#
# Design
# ------
# `plot_composition_scatter()` turns the user's own composition scatter
# (proportion of each cell state in WT on x, in the mutant on y, one panel
# per stage, y = x reference line, bold repelled labels) into a function,
# and replaces `plot_meta_scatter()`.
#
# * What is plotted is chosen by the arguments: each point is a level of
#   `group.by` (e.g. cell state); `split.by` is the column whose two levels
#   (`x`, `y`) give the axes (e.g. genotype); `facet.by` gives one panel per
#   level (e.g. stage). Swapping `split.by` and `facet.by` compares stages
#   within each genotype instead. (Here `split.by` picks the axes; it is not
#   a facet as in the other plots, where `facet.by` plays that role.)
# * Proportions are computed within each panel and condition: the share of
#   the cells of that condition (and panel) in each group, over all groups,
#   also those not shown (`groups`). A subset of cells (`cells`, e.g. one
#   lineage) changes the denominator to that subset.
# * Replicates (optional, `sample.by`): pooling cells gives one proportion
#   per condition, in which a large sample counts more than a small one.
#   With `sample.by`, a proportion is computed per sample and the point is
#   the mean of the samples of each condition, with error bars on both axes:
#   the horizontal bar from the samples of the x condition, the vertical one
#   from those of the y condition. `error` (2026-10-14): SEM (default), SD,
#   or a 95% confidence interval of the mean (t-based, mean +- t(0.975,
#   n - 1) * SEM, since conditions have few samples). Bars have caps of a
#   fixed physical size (as in GraphPad Prism) and are 0.75 pt black; an end
#   past the axis (e.g. below 0, or 0 on a log axis) is cut at the axis.
#   Without `sample.by` (default), the pooled proportion, as in the user's
#   code.
# * No statistics: whether a composition change is significant needs a test
#   on samples (e.g. propeller, scCODA, or Milo for neighbourhoods), not on
#   cells.
# * Axes: the same range on x and y and in every panel, square panels (F2),
#   so that the y = x line is the diagonal; `log = TRUE` spreads rare
#   populations, which are crushed near the origin on linear axes.


#' @title Count and proportion of cells per group and condition
#'
#' @description Counts the cells of each group (e.g. cell type) in each
#'              condition, and the proportion they make up of that
#'              condition's cells, optionally within panels (`facet.by`) and
#'              per sample (`sample.by`). The data behind
#'              \code{\link{plot_composition_scatter}}.
#'
#' @param object A Seurat object, or a data frame with one row per cell (e.g.
#'               its metadata)
#'
#' @param group.by Column of the groups counted (e.g. cell type)
#'
#' @param split.by Column of the conditions (e.g. genotype)
#'
#' @param facet.by Column within whose levels proportions are computed
#'                 separately (e.g. stage), or NULL (default)
#'
#' @param sample.by Column of the samples, to count each sample separately,
#'                  or NULL (default) to pool cells
#'
#' @param cells Cells (rows) to use, as names or indices. If NULL (default),
#'              all
#'
#' @return A data frame with one row per group and condition (and facet and
#'         sample): `facet`, `sample` (when given), `condition`, `group`,
#'         `count`, `total` (cells of the condition) and `prop`. Groups with
#'         no cells in a condition are included with a count of 0. Cells
#'         without a value in any of the columns are left out, with a message
#'
#' @examples
#' meta <- data.frame(
#'     cell_type = c("A", "A", "B", "B", "B", "C"),
#'     genotype = c("WT", "Mut", "WT", "Mut", "Mut", "WT")
#' )
#' summarize_group_meta(meta, group.by = "cell_type", split.by = "genotype")
#'
#' @seealso \code{\link{plot_composition_scatter}}
#'
#' @export
#'
summarize_group_meta <- function(object,
                                 group.by,
                                 split.by,
                                 facet.by = NULL,
                                 sample.by = NULL,
                                 cells = NULL) {
    if (!(is.data.frame(object) || inherits(object, "Seurat"))) {
        stop("`object` must be a Seurat object or a data frame", call. = FALSE)
    }
    vars <- c(group = group.by, condition = split.by, facet = facet.by, sample = sample.by)
    for (v in c("group.by", "split.by", "facet.by", "sample.by")) {
        val <- get(v)
        if (!is.null(val) && !(is.character(val) && length(val) == 1)) {
            stop("`", v, "` must be a single column name", call. = FALSE)
        }
    }
    if (is.data.frame(object) && is.null(rownames(object))) rownames(object) <- seq_len(nrow(object))
    cells <- plot_cells(object, cells)
    df <- as.data.frame(lapply(vars, function(v) fetch_plot_group(object, v, cells)))
    names(df) <- names(vars)
    keep <- stats::complete.cases(df)
    if (!all(keep)) {
        message(sum(!keep), " cell(s) without a value of `", paste(vars, collapse = "`, `"), "` are left out")
        df <- df[keep, , drop = FALSE]
    }
    if (nrow(df) == 0) stop("No cells to count", call. = FALSE)
    df[] <- lapply(df, droplevels)

    # Counts for every group in every condition (and facet / sample) present
    strata <- setdiff(names(df), "group")
    present <- unique(df[strata])
    grid <- merge(present, data.frame(group = factor(levels(df$group), levels = levels(df$group))))
    counts <- stats::aggregate(list(count = rep(1L, nrow(df))), df, length)
    out <- merge(grid, counts, all.x = TRUE)
    out$count[is.na(out$count)] <- 0L
    totals <- stats::aggregate(list(total = out$count), out[strata], sum)
    out <- merge(out, totals)
    out$prop <- out$count / out$total
    ord <- do.call(order, out[c(intersect(c("facet", "sample", "condition"), names(out)), "group")])
    out <- out[ord, c(intersect(c("facet", "sample", "condition"), names(out)), "group", "count", "total", "prop")]
    rownames(out) <- NULL
    out
}


#' @title Composition scatter: group proportions in two conditions
#'
#' @description Plots, for each group of cells (e.g. cell type or state), its
#'              proportion (or number of cells) in one condition on x against
#'              another condition on y, with the y = x line: groups above the
#'              line are more abundant in the y condition. Optionally one
#'              panel per level of `facet.by` (e.g. stage), and per-sample
#'              proportions with their mean and SEM (`sample.by`). Square
#'              panels of a fixed size; save it with \code{\link{save_plot}};
#'              the table behind it is returned by
#'              \code{\link{export_source_data}}.
#'
#' @param object A Seurat object, or a data frame with one row per cell (e.g.
#'               its metadata)
#'
#' @param group.by Column of the groups, one point each (e.g. cell state)
#'
#' @param split.by Column of the conditions; its levels `x` and `y` are
#'                 compared on the axes (e.g. genotype)
#'
#' @param x,y Levels of `split.by` on the x- and y-axis. If NULL (default),
#'            the first and second level
#'
#' @param facet.by Column giving one panel per level (e.g. stage), with
#'                 proportions computed within each. Default is NULL
#'
#' @param sample.by Column of the samples. If given, proportions are computed
#'                  per sample and each point is the mean of the samples of a
#'                  condition, with error bars on both axes (see `error`), so
#'                  that every sample counts the same. If NULL (default),
#'                  cells of all samples are pooled into one proportion per
#'                  condition
#'
#' @param error Error bars with `sample.by`: "sem" (default; standard error
#'              of the mean), "sd" (standard deviation of the samples) or
#'              "ci95" (95\% confidence interval of the mean, mean +-
#'              t(0.975, n - 1) x SEM). The horizontal bar is computed from
#'              the samples of the x condition, the vertical one from those
#'              of the y condition. A condition with one sample has no bar
#'
#' @param cells Cells (rows) to use, as names or indices, e.g. one lineage:
#'              proportions are then of that subset. If NULL (default), all
#'
#' @param groups Groups to draw. Proportions are still computed over all
#'               groups. If NULL (default), all
#'
#' @param value "percent" (default; the proportion, labelled in \%) or
#'              "count" (the number of cells)
#'
#' @param log Whether to draw both axes on a log10 scale, default is FALSE.
#'            Spreads rare groups; groups with 0 cells in a condition cannot
#'            be drawn and are left out, with a message
#'
#' @param fold.lines Fold changes drawn as dashed lines around y = x, e.g. 2
#'                   draws y = 2x and y = x / 2. Default is NULL (none)
#'
#' @param label Which points to name (needs ggrepel): TRUE (default) for
#'              all, FALSE for none, or the groups to name. Labels of many
#'              rare groups crowd near the origin on linear axes; name a few,
#'              or use `log = TRUE`
#'
#' @param label.size Size of the labels in mm. If NULL (default), the small
#'                   text size of the theme (9 pt)
#'
#' @param cols Colours of the groups: NULL (default) for the colours
#'             \code{\link{DimPlot_bioToolkit}} gives them, a qualitative
#'             palette name, or colours named by group
#'
#' @param pt.size Size of the points, default is 4
#'
#' @param panel.width,panel.height Size of each (square) panel in inches,
#'                                 default 2.5; `panel.height` defaults to
#'                                 `panel.width`
#'
#' @param legend.position Position of the group legend, default is "none"
#'                        (points are named)
#'
#' @param bg Background colour, default is "transparent"
#'
#' @details
#' The source data (\code{\link{export_source_data}}) hold one row per group
#' (and panel): the values on x and y, the cell counts and totals of both
#' conditions and, with `sample.by`, the SEM, SD, number of samples and the
#' half-width of the error bars drawn (`x_error`, `y_error`).
#' Per-sample proportions are returned by \code{\link{summarize_group_meta}}.
#'
#' Groups with no cells in either condition (of a panel, or among `cells`)
#' are not drawn.
#'
#' The plot shows where proportions differ, not whether the difference is
#' significant: that needs a test on samples, not cells (e.g. propeller,
#' scCODA).
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' # Pooled proportions of each cell state, WT vs Mut, one panel per stage
#' plot_composition_scatter(obj, group.by = "cell_state", split.by = "genotype",
#'                          x = "WT", y = "Mut", facet.by = "stage")
#'
#' # Per-sample proportions: mean +- SEM; log axes for rare states
#' p <- plot_composition_scatter(obj, group.by = "cell_state", split.by = "genotype",
#'                               x = "WT", y = "Mut", sample.by = "sample",
#'                               log = TRUE, fold.lines = 2)
#' save_plot(p, "composition.pdf")
#' }
#'
#' @seealso \code{\link{summarize_group_meta}}, \code{\link{save_plot}}
#'
#' @import ggplot2
#' @export
#'
plot_composition_scatter <- function(
    object,
    group.by,
    split.by,
    x = NULL,
    y = NULL,
    facet.by = NULL,
    sample.by = NULL,
    error = c("sem", "sd", "ci95"),
    cells = NULL,
    groups = NULL,
    value = c("percent", "count"),
    log = FALSE,
    fold.lines = NULL,
    label = TRUE,
    label.size = NULL,
    cols = NULL,
    pt.size = 4,
    panel.width = 2.5,
    panel.height = panel.width,
    legend.position = "none",
    bg = "transparent"
) {
    value <- match.arg(value)
    error <- match.arg(error)
    if (!(isTRUE(label) || isFALSE(label) || is.character(label))) {
        stop("`label` must be TRUE, FALSE or the groups to name", call. = FALSE)
    }
    label_on <- !isFALSE(label) && length(label) > 0
    if (label_on) check_installed("ggrepel", reason = "to label the points")
    sz <- theme_sizes()
    if (is.null(label.size)) label.size <- pt_to_size(sz$small)

    tab <- summarize_group_meta(object, group.by, split.by, facet.by, sample.by, cells)
    conditions <- levels(tab$condition)
    if (length(conditions) < 2) stop("`split.by` has fewer than two levels in the cells used", call. = FALSE)
    if (is.null(x)) x <- conditions[1]
    if (is.null(y)) y <- setdiff(conditions, x)[1]
    for (lv in c(x, y)) {
        if (!(lv %in% conditions)) {
            stop("'", lv, "' is not a level of `", split.by, "`: ", paste(conditions, collapse = ", "), call. = FALSE)
        }
    }
    if (x == y) stop("`x` and `y` must be different levels", call. = FALSE)
    message("x: ", x, ", y: ", y, " (levels of `", split.by, "`)")

    wide <- composition_wide(tab, x, y, value)
    # A group with no cells in either condition (e.g. a myeloid cell type in a
    # lymphoid panel of `facet.by`, or outside `cells`) is not present there:
    # it is not drawn, rather than drawn and labelled at the origin
    wide <- wide[wide$x_value > 0 | wide$y_value > 0, , drop = FALSE]
    # Groups drawn; colours from all levels in the object, as DimPlot
    group_all <- fetch_plot_group(object, group.by, plot_cells(object, NULL))
    if (!is.null(groups)) {
        missing <- setdiff(groups, levels(wide$group))
        if (length(missing)) warning("`groups` not found: ", paste(missing, collapse = ", "), call. = FALSE)
        wide <- wide[wide$group %in% groups, , drop = FALSE]
    }
    if (log) {
        zero <- wide$x_value <= 0 | wide$y_value <= 0
        if (any(zero)) {
            message(sum(zero), " group(s) with no cells in a condition cannot be drawn on log axes: ",
                    paste(unique(wide$group[zero]), collapse = ", "))
            wide <- wide[!zero, , drop = FALSE]
        }
    }
    if (nrow(wide) == 0) stop("No groups to draw", call. = FALSE)
    group_cols <- group_colours(cols, group_all)$cols

    # Half-width of the error bars (per sample): SEM, SD or t-based 95% CI,
    # each axis from the samples of its own condition
    per_sample <- !is.null(wide$x_sem)
    if (per_sample) {
        for (a in c("x", "y")) {
            n <- wide[[paste0(a, "_n_samples")]]
            sem <- wide[[paste0(a, "_sem")]]
            wide[[paste0(a, "_error")]] <- switch(
                error,
                sem = sem,
                sd = wide[[paste0(a, "_sd")]],
                ci95 = ifelse(n > 1, stats::qt(0.975, pmax(n - 1, 1)) * sem, 0)
            )
        }
    }

    # Same range on both axes and in every panel
    ends <- c(wide$x_value, wide$y_value)
    if (per_sample) ends <- c(ends, wide$x_value + wide$x_error, wide$y_value + wide$y_error)
    lim <- if (log) {
        r <- range(ends[ends > 0])
        c(r[1] / 1.25, r[2] * 1.25)
    } else {
        c(0, max(ends) * 1.05)
    }
    labels <- if (value == "percent") scales::label_percent(drop0trailing = TRUE) else
        scales::label_number(scale_cut = scales::cut_short_scale(), drop0trailing = TRUE)
    axis_scale <- function(fun_lin, fun_log) {
        if (log) fun_log(limits = lim, labels = labels, expand = expansion(0)) else
            fun_lin(limits = lim, labels = labels, expand = expansion(0))
    }
    what <- if (value == "percent") "% of cells" else "Cells"

    lw <- pt_to_linewidth(sz$line)
    p <- ggplot(wide, aes(x = .data$x_value, y = .data$y_value)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey50",
                    linewidth = pt_to_linewidth(0.75))
    for (k in fold.lines) {
        if (!(is.numeric(k) && k > 0 && k != 1)) stop("`fold.lines` must be positive numbers other than 1",
                                                     call. = FALSE)
        # On log axes the lines y = kx are parallel to the diagonal
        slope_lines <- if (log) {
            list(geom_abline(slope = 1, intercept = log10(k), linetype = "dashed", colour = "grey70",
                             linewidth = lw),
                 geom_abline(slope = 1, intercept = -log10(k), linetype = "dashed", colour = "grey70",
                             linewidth = lw))
        } else {
            list(geom_abline(slope = k, intercept = 0, linetype = "dashed", colour = "grey70", linewidth = lw),
                 geom_abline(slope = 1 / k, intercept = 0, linetype = "dashed", colour = "grey70", linewidth = lw))
        }
        p <- p + slope_lines
    }
    if (per_sample) {
        # Bar ends are kept within the axis (an SD below 0, or 0 on a log
        # axis, is cut there). Caps are 0.1 in wide whatever the axis range:
        # `width` is in the units of the axis drawn (log10 units on a log
        # axis), with no expansion
        bars <- wide[(wide$x_error > 0) | (wide$y_error > 0), , drop = FALSE]
        bars$x_low <- pmax(bars$x_value - bars$x_error, lim[1])
        bars$x_high <- bars$x_value + bars$x_error
        bars$y_low <- pmax(bars$y_value - bars$y_error, lim[1])
        bars$y_high <- bars$y_value + bars$y_error
        span <- if (log) diff(log10(lim)) else diff(lim)
        cap <- 0.1 * span
        bar_lw <- pt_to_linewidth(0.75)
        p <- p +
            geom_errorbar(
                aes(xmin = .data$x_low, xmax = .data$x_high),
                data = bars[bars$x_error > 0, , drop = FALSE],
                width = cap / panel.height,
                colour = "black",
                linewidth = bar_lw,
                orientation = "y"
            ) +
            geom_errorbar(
                aes(ymin = .data$y_low, ymax = .data$y_high),
                data = bars[bars$y_error > 0, , drop = FALSE],
                width = cap / panel.width,
                colour = "black",
                linewidth = bar_lw
            )
    }
    p <- p +
        geom_point(aes(fill = .data$group), shape = 21, colour = "black", size = pt.size,
                   stroke = pt_to_stroke(0.75)) +
        scale_fill_manual(values = group_cols, drop = TRUE, name = group.by) +
        axis_scale(scale_x_continuous, scale_x_log10) +
        axis_scale(scale_y_continuous, scale_y_log10) +
        labs(x = paste(what, "in", x), y = paste(what, "in", y))
    if (label_on) {
        label_groups <- if (isTRUE(label)) levels(wide$group) else label
        p <- p + ggrepel::geom_text_repel(
            aes(label = .data$group), size = label.size, fontface = "bold", colour = "black",
            data = wide[wide$group %in% label_groups, , drop = FALSE],
            box.padding = 0.5, point.padding = 0.4, min.segment.length = 0, segment.colour = "grey50",
            segment.size = lw, max.overlaps = Inf, seed = 42, max.time = 2
        )
    }
    n_panels <- 1L
    if (!is.null(facet.by)) {
        p <- p + facet_wrap(vars(.data$facet), nrow = 1, axes = "all")
        n_panels <- length(unique(wide$facet))
    }
    p <- p +
        theme_bioToolkit("axis", bg = bg) +
        panel_size_theme(panel.width, panel.height, ncol = n_panels) +
        coord_cartesian(clip = "off") +
        theme(legend.position = legend.position,
              strip.text = element_text(size = sz$large, face = "bold"))
    attr(p, "source_data") <- wide
    p
}


#' One row per group with its values in the two conditions
#'
#' @param tab Result of \code{summarize_group_meta()}
#'
#' @param x,y Conditions on the axes
#'
#' @param value "percent" or "count"
#'
#' @return A data frame: `facet` (if any), `group`, `x_value`, `y_value`, and
#'         `x_count`, `x_total`, `y_count`, `y_total` (pooled), or `x_sd`,
#'         `x_sem`, `y_sd`, `y_sem`, `x_n_samples`, `y_n_samples` (per
#'         sample)
#'
#' @keywords internal
composition_wide <- function(tab,
                             x,
                             y,
                             value) {
    v <- if (value == "percent") tab$prop else tab$count
    tab$v <- v
    keys <- intersect("facet", names(tab))
    per_sample <- "sample" %in% names(tab)
    side <- function(cond, prefix) {
        d <- tab[tab$condition == cond, , drop = FALSE]
        if (per_sample) {
            # Mean, SD and SEM over the samples of the condition
            agg <- stats::aggregate(list(value = d$v), d[c(keys, "group")], mean)
            agg$sd <- stats::aggregate(list(s = d$v), d[c(keys, "group")],
                                       function(z) if (length(z) > 1) stats::sd(z) else 0)$s
            agg$sem <- stats::aggregate(list(s = d$v), d[c(keys, "group")],
                                        function(z) if (length(z) > 1) stats::sd(z) / sqrt(length(z)) else 0)$s
            agg$n_samples <- stats::aggregate(list(n = d$v), d[c(keys, "group")], length)$n
        } else {
            agg <- d[c(keys, "group")]
            agg$value <- d$v
            agg$count <- d$count
            agg$total <- d$total
        }
        names(agg)[!(names(agg) %in% c(keys, "group"))] <- paste0(prefix, "_",
                                                                  names(agg)[!(names(agg) %in% c(keys, "group"))])
        agg
    }
    out <- merge(side(x, "x"), side(y, "y"), by = c(keys, "group"))
    if (nrow(out) == 0) stop("No group has cells in both '", x, "' and '", y, "'", call. = FALSE)
    out$group <- factor(as.character(out$group), levels = levels(tab$group))
    if ("facet" %in% keys) out$facet <- factor(as.character(out$facet), levels = levels(tab$facet))
    out <- out[do.call(order, out[c(keys, "group")]), , drop = FALSE]
    rownames(out) <- NULL
    out
}
