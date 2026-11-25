# Cell composition between two conditions


#' @title Count and proportion of cells per group and condition
#'
#' @description Counts the cells (or genes, peaks...) of each group (e.g. cell
#'              type) in each condition, and the proportion they make up of
#'              that condition's cells, optionally within panels (`facet.by`)
#'              and per sample (`sample.by`). The data behind
#'              \code{\link{plot_composition_scatter}} and
#'              \code{\link{plot_composition_bars}}.
#'
#' @param object A Seurat object, or a table with one row per cell, gene,
#'               peak... (its metadata, a gene or peak annotation): anything
#'               \code{as.data.frame()} converts, e.g. a data frame, tibble,
#'               data.table, Bioconductor \code{DataFrame} or \code{GRanges}
#'
#' @param group.by Column of the groups counted (e.g. cell type)
#'
#' @param split.by Column of the conditions (e.g. genotype), or NULL (default)
#'                 to count all cells as one condition, "All <unit>"
#'
#' @param facet.by Column within whose levels proportions are computed
#'                 separately (e.g. stage), or NULL (default)
#'
#' @param sample.by Column of the samples, to count each sample separately,
#'                  or NULL (default) to pool cells
#'
#' @param cells Cells (rows) to use, as names or indices. If NULL (default),
#'              all. A table without row names of its own is numbered, and
#'              then takes row indices only
#'
#' @param unit What a row is, used in messages and in the name of the single
#'             condition when `split.by` is NULL: "cells" (default), "genes",
#'             "peaks"...
#'
#' @details Columns must be categorical: factors, characters, logicals, or
#' whole numbers (e.g. cluster numbers), which become factors in increasing
#' order. Other numbers must be binned first, e.g. with \code{cut()}.
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
#' # One condition: the share of all rows in each group
#' summarize_group_meta(meta, group.by = "cell_type", unit = "genes")
#'
#' @seealso \code{\link{plot_composition_scatter}},
#'          \code{\link{plot_composition_bars}}
#'
#' @export
#'
summarize_group_meta <- function(object,
                                 group.by,
                                 split.by = NULL,
                                 facet.by = NULL,
                                 sample.by = NULL,
                                 cells = NULL,
                                 unit = "cells") {
    object <- as_plot_table(object)
    if (!(is.character(unit) && length(unit) == 1)) {
        stop("`unit` must be a single word, e.g. \"cells\"", call. = FALSE)
    }
    for (v in c("group.by", "split.by", "facet.by", "sample.by")) {
        val <- get(v)
        if (!is.null(val) && !(is.character(val) && length(val) == 1)) {
            stop("`", v, "` must be a single column name", call. = FALSE)
        }
    }
    vars <- c(
        group = group.by,
        condition = split.by,
        facet = facet.by,
        sample = sample.by
    )

    # Categories only: a column of measurements would give one group per value
    meta <- if (is.data.frame(object)) object else object[[]]
    for (v in names(vars)) {
        if (vars[[v]] %in% colnames(meta)) check_categorical(meta[[vars[[v]]]], vars[[v]])
    }

    # Numbered rows have no names to select them by
    if (isTRUE(attr(object, "row_numbers")) && is.character(cells)) {
        stop(
            "The table has no row names of its own (its rows are numbered); ",
            "give `cells` as row indices",
            call. = FALSE
        )
    }
    cells <- plot_cells(object, cells)
    df <- as.data.frame(lapply(vars, function(v) fetch_plot_group(object, v, cells)))
    names(df) <- names(vars)
    # No split.by: every row in one condition. (fetch_plot_group() would give
    # the identities of a Seurat object.)
    if (is.null(split.by)) {
        df$condition <- factor(rep(paste("All", unit), length(cells)))
        df <- df[c("group", "condition", intersect(c("facet", "sample"), names(df)))]
    }
    keep <- stats::complete.cases(df)
    if (!all(keep)) {
        columns <- paste(vars, collapse = "`, `")
        message(sum(!keep), " ", unit, " without a value of `", columns, "` are left out")
        df <- df[keep, , drop = FALSE]
    }
    if (nrow(df) == 0) stop("No ", unit, " to count", call. = FALSE)
    df[] <- lapply(df, droplevels)

    # Counts for every group in every condition (and facet / sample) present
    strata <- setdiff(names(df), "group")
    present <- unique(df[strata])
    all_groups <- factor(levels(df$group), levels = levels(df$group))
    grid <- merge(present, data.frame(group = all_groups))
    counts <- stats::aggregate(list(count = rep(1L, nrow(df))), df, length)
    out <- merge(grid, counts, all.x = TRUE)
    out$count[is.na(out$count)] <- 0L
    totals <- stats::aggregate(list(total = out$count), out[strata], sum)
    out <- merge(out, totals)
    out$prop <- out$count / out$total
    keys <- c(intersect(c("facet", "sample", "condition"), names(out)), "group")
    ord <- do.call(order, out[keys])
    out <- out[ord, c(keys, "count", "total", "prop")]
    rownames(out) <- NULL
    out
}


#' A Seurat object, or a table as a data frame
#'
#' Anything that is not a Seurat object is converted with
#' \code{as.data.frame()} (data frames, tibbles, data.tables, Bioconductor
#' \code{DataFrame}s, \code{GRanges}). Tables without row names of their own
#' (automatic row names, or a \code{GRanges} with duplicated names) are
#' numbered, and flagged with the attribute `row_numbers`.
#'
#' @param object A Seurat object or a table
#'
#' @return The Seurat object, or a data frame with unique row names
#'
#' @keywords internal
as_plot_table <- function(object) {
    if (inherits(object, "Seurat")) return(object)
    # Already converted (e.g. by a plot function before summarize_group_meta()):
    # converting again would lose the `row_numbers` flag
    if (is.data.frame(object) && !is.null(attr(object, "row_numbers"))) return(object)
    convert <- function(x) tryCatch(as.data.frame(x), error = function(e) NULL)
    out <- convert(object)
    numbered <- FALSE
    # A GRanges with duplicated range names cannot become a data frame with
    # them as row names: drop the names
    if (is.null(out) && methods::is(object, "GenomicRanges")) {
        names(object) <- NULL
        out <- convert(object)
        numbered <- TRUE
    }
    if (!is.data.frame(out)) {
        stop(
            "`object` must be a Seurat object or a table that as.data.frame() can convert",
            call. = FALSE
        )
    }
    # Automatic row names (tibbles, data frames read from files) are numbers
    if (.row_names_info(out) < 0) {
        rownames(out) <- as.character(seq_len(nrow(out)))
        numbered <- TRUE
    }
    attr(out, "row_numbers") <- numbered
    out
}


#' Stop if a column is not categorical
#'
#' Factors, characters and logicals are categorical; so are whole numbers
#' (cluster or chromosome numbers). Other numbers are measurements.
#'
#' @param values The column's values
#'
#' @param column The column's name, for the message
#'
#' @return `values`, invisibly
#'
#' @keywords internal
check_categorical <- function(values,
                              column) {
    if (is.numeric(values)) {
        v <- values[is.finite(values)]
        if (any(v != round(v))) {
            stop(
                "`", column, "` is not categorical (it has non-integer numbers); ",
                "bin it with cut() first",
                call. = FALSE
            )
        }
    }
    invisible(values)
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
#' @param object A Seurat object, or a table with one row per cell (e.g. its
#'               metadata): anything \code{as.data.frame()} converts
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
    # Seurat object or any table (as.data.frame())
    object <- as_plot_table(object)
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


#' @title Stacked composition bars
#'
#' @description Stacked bars of how the cells (or genes, peaks...) of each
#'              bar are split among groups: one bar per level of `split.by`
#'              (e.g. condition or sample), one coloured segment per level of
#'              `group.by` (e.g. cell type), optionally one panel per level of
#'              `facet.by`. Neighbouring bars of a panel are joined by bands
#'              in a lighter tone of each group's colour (`connect`), so that a
#'              shift in a group reads at a glance. Bars have a fixed width, gap and
#'              length in inches, so that a figure keeps its proportions
#'              whatever the number of bars; save it with
#'              \code{\link{save_plot}}; the table behind it is returned by
#'              \code{\link{export_source_data}}.
#'
#' @param object A Seurat object, or a table with one row per cell, gene,
#'               peak... with categorical columns (e.g. its metadata, a gene
#'               or peak annotation): anything \code{as.data.frame()}
#'               converts, e.g. a data frame, tibble, data.table,
#'               Bioconductor \code{DataFrame} or \code{GRanges}
#'
#' @param group.by Column of the groups, the coloured segments (e.g. cell
#'                 type)
#'
#' @param split.by Column giving one bar per level (e.g. condition or
#'                 sample), or NULL (default) for one bar of all cells. Each
#'                 bar's segments are shares of that bar's cells
#'
#' @param facet.by Column giving one panel per level (e.g. condition when the
#'                 bars are samples), or NULL (default). A panel holds only
#'                 the bars with cells in it
#'
#' @param value What the value axis shows: "percent" (default; 0-100\%),
#'              "fraction" (0-1) or "count" (the number of cells; bars of
#'              different lengths)
#'
#' @param unit What a row is, for the axis title and messages: "cells"
#'             (default), "genes", "peaks"...
#'
#' @param cells Cells (rows) to use, as names or indices, e.g. one lineage:
#'              shares are then of that subset. If NULL (default), all
#'
#' @param order.by A group: the bars of each panel are sorted by decreasing
#'                 share of that group. If NULL (default), bars follow the
#'                 levels of `split.by`
#'
#' @param flip Whether to draw horizontal bars (first bar on top). Default is
#'             FALSE (vertical bars)
#'
#' @param connect Whether and how to join neighbouring bars of a panel with
#'                bands in a light tone of each group's colour: "none"
#'                (default; no bands), "curve" (smooth, as in alluvial plots)
#'                or "straight"
#'
#' @param connect.lighten How light the bands are: the share of white mixed
#'                        into each group's colour, from 0 (the group's
#'                        colour) to 1 (white). Default 0.8, light enough to
#'                        stay behind the bars
#'
#' @param label Whether to write each segment's value inside it (black, or
#'              white on a dark segment), where the text fits in the
#'              segment. Default is FALSE
#'
#' @param label.min Smallest share of a bar whose segment is labelled,
#'                  default 0.05
#'
#' @param show.n Whether to write the number of cells of each bar beyond its
#'               end. Default is FALSE
#'
#' @param value.axis.line Whether to draw the line of the value axis (\%,
#'                        fraction or count). Default is FALSE: its ticks and
#'                        labels only
#'
#' @param cols Colours of the groups: NULL (default) for the colours
#'             \code{\link{DimPlot_bioToolkit}} gives them, a qualitative
#'             palette name, or colours named by group
#'
#' @param facet.cols Colours of the bars under the panel names: NULL
#'                   (default) for black, a qualitative palette name, or
#'                   colours named by level of `facet.by`
#'
#' @param bar.width Width of a bar in inches, default 0.25
#'
#' @param bar.gap Gap between neighbouring bars in inches, where the bands
#'                are drawn; default 0.125
#'
#' @param bar.length Length of the value axis in inches (a whole bar in
#'                   percent or fraction), default 2.5
#'
#' @param border.width Width of the black borders around the bars and
#'                     between their segments, in points; default 1. 0 for
#'                     no borders
#'
#' @param panel.spacing Gap between panels in inches, default 0.2
#'
#' @param legend.position Position of the legend, one of "right" (default),
#'                        "left", "top", "bottom", "inside" or "none"
#'
#' @param legend.ncol Number of legend columns. If NULL (default), ggplot2
#'                    chooses based on the number of groups and the position
#'
#' @param legend.title Title of the legend. If NULL (default), the name of
#'                     `group.by`; "" for no title
#'
#' @param bg Background colour, default is "transparent"
#'
#' @details
#' Shares are within each bar: the segments of a bar add up to 100\%. To see
#' instead which bars make up each group (e.g. which samples make up each
#' cluster, a check of integration), swap `group.by` and `split.by`.
#'
#' Pooling the cells of several samples into one bar lets a large sample
#' count more than a small one; `split.by = "sample"` with
#' `facet.by = "condition"` shows every sample.
#'
#' The plot shows where composition differs, not whether the difference is
#' significant: that needs a test on samples, not cells (e.g. propeller,
#' scCODA).
#'
#' The source data (\code{\link{export_source_data}}) hold one row per
#' segment: the panel, bar, group, its `count`, the bar's `total`, the share
#' `prop`, the `value` drawn and the segment's ends on the value axis
#' (`lower`, `upper`).
#'
#' @return A ggplot object
#'
#' @examples
#' meta <- data.frame(
#'     cell_type = rep(c("T", "B", "NK", "T", "B", "Mono"), times = c(30, 10, 5, 20, 15, 10)),
#'     sample = rep(c("s1", "s2"), times = c(45, 45)),
#'     condition = rep(c("Control", "Treated"), times = c(45, 45))
#' )
#' plot_composition_bars(meta, group.by = "cell_type", split.by = "sample")
#'
#' \dontrun{
#' # One bar per sample, one panel per condition
#' p <- plot_composition_bars(
#'     obj,
#'     group.by = "cell_type",
#'     split.by = "sample",
#'     facet.by = "condition",
#'     show.n = TRUE
#' )
#' save_plot(p, "composition_bars.pdf")
#'
#' # Any table: the direction of genes per cell type
#' plot_composition_bars(de, group.by = "direction", split.by = "cell_type", unit = "genes")
#' }
#'
#' @seealso \code{\link{summarize_group_meta}},
#'          \code{\link{plot_composition_scatter}}, \code{\link{save_plot}}
#'
#' @import ggplot2
#' @export
#'
plot_composition_bars <- function(
    object,
    group.by,
    split.by = NULL,
    facet.by = NULL,
    value = c("percent", "fraction", "count"),
    unit = "cells",
    cells = NULL,
    order.by = NULL,
    flip = FALSE,
    connect = c("none", "curve", "straight"),
    connect.lighten = 0.8,
    label = FALSE,
    label.min = 0.05,
    show.n = FALSE,
    value.axis.line = FALSE,
    cols = NULL,
    facet.cols = NULL,
    bar.width = 0.25,
    bar.gap = 0.125,
    bar.length = 2.5,
    border.width = 1,
    panel.spacing = 0.2,
    legend.position = c("right", "left", "top", "bottom", "inside", "none"),
    legend.ncol = NULL,
    legend.title = NULL,
    bg = "transparent"
) {
    value <- match.arg(value)
    connect <- match.arg(connect)
    legend.position <- rlang::arg_match(legend.position)
    ncol_ok <- is.numeric(legend.ncol) && length(legend.ncol) == 1 && legend.ncol >= 1
    if (!is.null(legend.ncol) && !ncol_ok) {
        stop("`legend.ncol` must be NULL or a positive whole number", call. = FALSE)
    }
    title_ok <- is.character(legend.title) && length(legend.title) == 1
    if (!is.null(legend.title) && !title_ok) {
        stop("`legend.title` must be NULL or a single string", call. = FALSE)
    }
    for (arg in c("flip", "label", "show.n")) {
        if (!(is.logical(get(arg)) && length(get(arg)) == 1 && !is.na(get(arg)))) {
            stop("`", arg, "` must be TRUE or FALSE", call. = FALSE)
        }
    }
    sizes <- list(
        bar.width = bar.width,
        bar.length = bar.length,
        bar.gap = bar.gap,
        panel.spacing = panel.spacing
    )
    for (arg in names(sizes)) {
        x <- sizes[[arg]]
        # Bars need a width and a length; gaps may be 0
        least_ok <- if (arg %in% c("bar.width", "bar.length")) x > 0 else x >= 0
        if (!(is.numeric(x) && length(x) == 1 && !is.na(x) && least_ok)) {
            kind <- if (arg %in% c("bar.width", "bar.length")) "positive" else "non-negative"
            stop("`", arg, "` must be a single ", kind, " number (inches)", call. = FALSE)
        }
    }
    border_ok <- is.numeric(border.width) && length(border.width) == 1
    if (!(border_ok && !is.na(border.width) && border.width >= 0)) {
        stop("`border.width` must be a single non-negative number (points)", call. = FALSE)
    }
    line_ok <- is.logical(value.axis.line) && length(value.axis.line) == 1
    if (!(line_ok && !is.na(value.axis.line))) {
        stop("`value.axis.line` must be TRUE or FALSE", call. = FALSE)
    }
    lighten_ok <- is.numeric(connect.lighten) && length(connect.lighten) == 1
    if (!(lighten_ok && !is.na(connect.lighten) && connect.lighten >= 0 && connect.lighten <= 1)) {
        stop("`connect.lighten` must be a number from 0 to 1", call. = FALSE)
    }
    min_ok <- is.numeric(label.min) && length(label.min) == 1 && label.min >= 0 && label.min <= 1
    if (!min_ok) {
        stop("`label.min` must be a number from 0 to 1", call. = FALSE)
    }
    if (connect != "none" && bar.gap == 0) {
        message("`bar.gap` is 0: no room for the bands between bars")
        connect <- "none"
    }
    sz <- theme_sizes()
    text_size <- pt_to_size(sz$small)

    # Counts and shares; colours from all levels in the object, as DimPlot,
    # so that a group keeps its colour in a subset
    object <- as_plot_table(object)
    tab <- summarize_group_meta(
        object,
        group.by,
        split.by,
        facet.by,
        cells = cells,
        unit = unit
    )
    group_all <- fetch_plot_group(object, group.by, plot_cells(object, NULL))
    group_cols <- group_colours(cols, group_all)$cols
    drawn <- levels(tab$group)
    if (length(drawn) > 30) {
        warning(
            "`", group.by, "` has ", length(drawn), " levels: ",
            "the colours and legend will be hard to read",
            call. = FALSE
        )
    }

    # Segments, in category (pos) / value coordinates. Horizontal bars: the
    # first bar on top, the first group at the left
    seg <- composition_segments(tab, value, order.by, first.at.start = flip)
    if (flip) seg$pos <- max(seg$pos) + 1 - seg$pos
    has_facet <- !is.null(facet.by)
    # A slot is a bar and its gap; one slot is one unit of the category axis
    slot <- bar.width + bar.gap
    half <- bar.width / slot / 2
    seg$cat_min <- seg$pos - half
    seg$cat_max <- seg$pos + half
    bars <- unique(seg[c(if (has_facet) "facet", "bar", "pos", "total")])

    # Width in inches of a label of `chars` characters at the small text size
    # (about 0.6 em per character; 0.65 in bold)
    label_in <- function(chars,
                         em = 0.6) {
        chars * em * sz$small / 72
    }

    # Each panel spans its bars and half a gap on either side, in slots. A
    # panel of few bars under a long name (vertical bars) is widened to the
    # name, keeping the bars' width: the extra room is shared on both sides.
    # Panels are then laid end to end on the category axis, so that a
    # widened panel never reaches the bars (and names) of the next one
    panel_of <- if (has_facet) bars$facet else factor(rep("all", nrow(bars)))
    panels <- split(seq_len(nrow(bars)), panel_of, drop = TRUE)
    panel_slots <- vapply(panels, function(rows) {
        slots <- length(rows)
        if (has_facet && !flip) {
            name <- as.character(bars$facet[rows[1]])
            name_in <- label_in(nchar(name), em = 0.65) + 4 / 72
            slots <- max(slots, name_in / slot)
        }
        slots
    }, numeric(1))
    # Panels in the order of their positions (reversed when horizontal)
    first_pos <- vapply(panels, function(rows) min(bars$pos[rows]), numeric(1))
    start <- 0.5
    extent <- list()
    for (k in order(first_pos)) {
        rows <- panels[[k]]
        # Shift the panel's bars to the middle of its slots
        first_bar <- start + (panel_slots[k] - length(rows)) / 2 + 0.5
        shift <- first_bar - min(bars$pos[rows])
        in_panel <- seg$pos %in% bars$pos[rows]
        seg$pos[in_panel] <- seg$pos[in_panel] + shift
        bars$pos[rows] <- bars$pos[rows] + shift
        e <- data.frame(cat = c(start, start + panel_slots[k]), val = 0)
        if (has_facet) e$facet <- bars$facet[rows[1]]
        extent[[k]] <- e
        start <- start + panel_slots[k]
    }
    # In the order of the panels (facet levels), as panel_size_theme() needs
    extent <- do.call(rbind, extent)
    panel_slots <- unname(panel_slots)
    seg$cat_min <- seg$pos - half
    seg$cat_max <- seg$pos + half

    # Value axis. With `show.n`, room beyond the longest bar for the numbers
    # (their height, or their width when they run along the axis), added to
    # the panel so that a whole bar stays `bar.length` long
    n_labels <- scales::label_comma()(bars$total)
    n_fits <- label_in(max(nchar(n_labels))) <= slot
    n_along <- flip || !n_fits
    room_in <- if (n_along) {
        label_in(max(nchar(n_labels))) + 4 / 72
    } else {
        1.4 * sz$small / 72
    }
    if (!show.n) room_in <- 0
    value_length <- bar.length + room_in
    unit_title <- paste0(toupper(substr(unit, 1, 1)), substring(unit, 2))
    value_title <- switch(
        value,
        percent = paste("% of", unit),
        fraction = paste("Fraction of", unit),
        count = unit_title
    )
    # Counts: the axis ends at the first break at or above the longest bar,
    # so that the axis line reaches past every bar
    # (whole numbers only: no "0.5 cells")
    longest <- max(seg$upper)
    count_breaks <- pretty(c(0, longest))
    count_breaks <- count_breaks[count_breaks == round(count_breaks)]
    if (max(count_breaks) < longest) count_breaks <- c(count_breaks, ceiling(longest))
    value_breaks <- if (value == "count") count_breaks else c(0, 0.25, 0.5, 0.75, 1)
    count_labels <- scales::label_number(
        scale_cut = scales::cut_short_scale(),
        drop0trailing = TRUE
    )
    value_labels <- switch(
        value,
        percent = scales::label_percent(),
        fraction = scales::label_number(drop0trailing = TRUE),
        count = count_labels
    )
    value_limits <- if (value == "count") range(count_breaks) else c(0, 1)
    # Room for the numbers in data units: `bar.length` inches per axis range
    room_add <- room_in / bar.length * diff(value_limits)
    value_scale <- function(fun) {
        fun(
            limits = value_limits,
            breaks = value_breaks,
            labels = value_labels,
            expand = expansion(add = c(0, room_add)),
            # Keep the cell numbers drawn just beyond the limit
            oob = scales::oob_keep,
            guide = guide_axis(cap = if (show.n) "upper" else "none")
        )
    }
    category_scale <- function(fun) {
        fun(
            breaks = bars$pos,
            labels = as.character(bars$bar),
            expand = expansion(0)
        )
    }

    # Aesthetics in category / value terms, mapped to x / y (swapped when
    # horizontal)
    cv <- function(cat, val) {
        if (flip) {
            aes(x = .data[[val]], y = .data[[cat]])
        } else {
            aes(x = .data[[cat]], y = .data[[val]])
        }
    }
    rect_aes <- if (flip) {
        aes(
            xmin = .data$lower,
            xmax = .data$upper,
            ymin = .data$cat_min,
            ymax = .data$cat_max,
            fill = .data$group
        )
    } else {
        aes(
            xmin = .data$cat_min,
            xmax = .data$cat_max,
            ymin = .data$lower,
            ymax = .data$upper,
            fill = .data$group
        )
    }

    p <- ggplot() +
        geom_blank(cv("cat", "val"), data = extent)
    if (connect != "none") {
        bands <- composition_bands(seg, half, shape = connect)
        if (nrow(bands) > 0) {
            # Solid lighter tones (not transparency), drawn under the bars
            band_cols <- unname(group_cols[as.character(bands$group)])
            bands$band_fill <- lighten_colour(band_cols, connect.lighten)
            band_aes <- cv("pos", "val")
            band_aes$group <- quote(.data$id)
            band_aes$fill <- quote(I(.data$band_fill))
            p <- p +
                geom_polygon(band_aes, data = bands, colour = NA, show.legend = FALSE)
        }
    }
    p <- p +
        geom_rect(
            rect_aes,
            data = seg,
            colour = if (border.width > 0) "black" else NA,
            linewidth = pt_to_linewidth(border.width)
        )

    # Segment values inside the segments: black, or white on a dark segment
    if (label) {
        lab <- seg[seg$prop >= label.min & seg$value > 0, , drop = FALSE]
        lab$mid <- (lab$lower + lab$upper) / 2
        lab$text <- switch(
            value,
            percent = paste0(round(lab$prop * 100), "%"),
            fraction = formatC(lab$prop, format = "f", digits = 2),
            count = as.character(lab$count)
        )
        # Only labels that fit inside their segment, with 2 pt to spare on
        # either end: along a horizontal bar the text's width, along a
        # vertical one its height (and its width across the bar)
        seg_in <- (lab$upper - lab$lower) / diff(value_limits) * bar.length
        text_w <- label_in(nchar(lab$text))
        text_h <- 1.2 * sz$small / 72
        pad <- 4 / 72
        fits <- if (flip) {
            seg_in >= text_w + pad & text_h <= bar.width
        } else {
            seg_in >= text_h + pad & text_w <= bar.width
        }
        lab <- lab[fits, , drop = FALSE]
        lab$text_colour <- text_contrast(group_cols[as.character(lab$group)])
        lab_aes <- cv("pos", "mid")
        lab_aes$label <- quote(.data$text)
        lab_aes$colour <- quote(I(.data$text_colour))
        p <- p +
            geom_text(lab_aes, data = lab, size = text_size, show.legend = FALSE)
    }

    # Number of cells beyond the end of each bar, in black
    if (show.n) {
        bars$end <- if (value == "count") bars$total else 1
        bars$n_label <- n_labels
        n_aes <- cv("pos", "end")
        n_aes$label <- quote(.data$n_label)
        n_text <- if (flip) {
            list(hjust = 0, vjust = 0.5, angle = 0)
        } else if (n_along) {
            list(hjust = 0, vjust = 0.5, angle = 90)
        } else {
            list(hjust = 0.5, vjust = 0, angle = 0)
        }
        p <- p +
            geom_text(
                n_aes,
                data = bars,
                size = text_size,
                colour = "black",
                hjust = n_text$hjust,
                vjust = n_text$vjust,
                angle = n_text$angle,
                nudge_x = if (flip) 0.02 * diff(value_limits) else 0,
                nudge_y = if (flip) 0 else 0.02 * diff(value_limits)
            )
    }

    # Panels: names in group bars, on the side away from the bar names
    if (has_facet) {
        facet_lv <- levels(droplevels(seg$facet))
        facet_fill <- if (is.null(facet.cols)) {
            stats::setNames(rep("black", length(facet_lv)), facet_lv)
        } else {
            facet_group <- factor(facet_lv, levels = facet_lv)
            group_colours(facet.cols, facet_group, arg = "facet.cols")$cols
        }
        # Each bar spans its whole panel (the layer adds half a unit at each end)
        panel_ends <- function(e) {
            data.frame(
                facet = e$facet[1],
                first = min(e$cat) + 0.5,
                last = max(e$cat) - 0.5
            )
        }
        ends <- do.call(rbind, lapply(split(extent, extent$facet, drop = TRUE), panel_ends))
        ends$fill <- unname(facet_fill[as.character(ends$facet)])
        side <- if (flip) "right" else "top"
        if (flip) {
            ends$y <- ends$first
            ends$yend <- ends$last
            p <- p +
                facet_grid(rows = vars(.data$facet), scales = "free_y", space = "free_y")
        } else {
            ends$x <- ends$first
            ends$xend <- ends$last
            p <- p +
                facet_grid(cols = vars(.data$facet), scales = "free_x", space = "free_x")
        }
        p <- p +
            group_bar_layer(ends, side = side, width = 3.5)
    }

    # Scales, titles and sizes
    if (flip) {
        p <- p +
            value_scale(scale_x_continuous) +
            category_scale(scale_y_continuous) +
            labs(x = value_title, y = split.by)
    } else {
        p <- p +
            category_scale(scale_x_continuous) +
            value_scale(scale_y_continuous) +
            labs(x = split.by, y = value_title)
    }
    # Legend: square keys outlined like the bars, the groups in stacking order
    fill_title <- if (is.null(legend.title)) group.by else legend.title
    p <- p +
        scale_fill_manual(values = group_cols, breaks = drawn, name = fill_title) +
        guides(fill = guide_legend(ncol = legend.ncol))
    panel_sizes <- if (flip) {
        panel_size_theme(
            value_length,
            panel_slots * slot,
            nrow = length(panel_slots),
            panel.spacing = panel.spacing
        )
    } else {
        panel_size_theme(
            panel_slots * slot,
            value_length,
            ncol = length(panel_slots),
            panel.spacing = panel.spacing
        )
    }
    # Bar names: level, or at 45 degrees when one does not fit its slot
    names_fit <- label_in(max(nchar(as.character(bars$bar)))) <= slot
    category_axis <- if (flip) {
        theme(axis.ticks.y = element_blank())
    } else if (names_fit) {
        theme(axis.ticks.x = element_blank())
    } else {
        theme(
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
        )
    }
    strip_text <- element_text(size = sz$small, face = "bold", colour = "black")
    p <- p +
        theme_bioToolkit("axis", bg = bg) +
        panel_sizes +
        category_axis +
        coord_cartesian(clip = "off") +
        theme(
            legend.position = legend.position,
            legend.key.size = unit(0.4, "cm"),
            legend.key.spacing.y = unit(0.08, "cm"),
            strip.placement = "outside",
            strip.clip = "off",
            strip.text.x.top = merge_element(
                element_text(margin = group_bar_margin("top", 3.5)),
                strip_text
            ),
            strip.text.y.right = merge_element(
                element_text(angle = 0, hjust = 0, margin = group_bar_margin("right", 3.5)),
                strip_text
            )
        )
    if (identical(legend.title, "")) {
        p <- p +
            theme(legend.title = element_blank())
    }
    # The value axis: ticks and labels, without its line unless asked
    if (!value.axis.line) {
        no_line <- if (flip) {
            theme(axis.line.x = element_blank())
        } else {
            theme(axis.line.y = element_blank())
        }
        p <- p +
            no_line
    }
    # Horizontal bars: the last value label ("100%") hangs past the panel's
    # end; keep a legend on the right clear of it
    if (flip && legend.position == "right") {
        p <- p +
            theme(legend.box.spacing = unit(0.25, "in"))
    }
    # A legend above or below the plot: its title above the keys
    if (legend.position %in% c("top", "bottom")) {
        p <- p +
            theme(legend.title.position = "top")
    }
    source <- seg[setdiff(names(seg), c("pos", "cat_min", "cat_max"))]
    rownames(source) <- NULL
    attr(p, "source_data") <- source
    p
}



#' Segments of stacked composition bars
#'
#' One row per segment (a group in a bar), placed and stacked. Bars are the
#' conditions of \code{summarize_group_meta()}; each panel (`facet`) holds
#' only the bars with cells in it. Positions run on over the whole plot (the
#' first panel's bars are 1..n, the next panel's n + 1...), so that a bar
#' name repeated in two panels gets its own place in each.
#'
#' @param tab Result of \code{summarize_group_meta()}
#'
#' @param value "percent" or "fraction" (the proportion) or "count"
#'
#' @param order.by A group: the bars of each panel are sorted by decreasing
#'                 proportion of that group (ties keep their order). NULL
#'                 keeps the order of the levels
#'
#' @param first.at.start Whether the first group is stacked at the start of
#'                       the value axis (left of a horizontal bar). FALSE puts
#'                       it at the end (top of a vertical bar), so that the
#'                       stack reads in the order of the legend
#'
#' @return A data frame: `facet` (if any), `bar`, `pos`, `group`, `count`,
#'         `total`, `prop`, `value` (drawn), `lower` and `upper` (ends of the
#'         segment on the value axis)
#'
#' @keywords internal
composition_segments <- function(tab,
                                 value = c("percent", "fraction", "count"),
                                 order.by = NULL,
                                 first.at.start = FALSE) {
    value <- match.arg(value)
    order_ok <- length(order.by) == 1 && order.by %in% levels(tab$group)
    if (!is.null(order.by) && !order_ok) {
        stop("`order.by` must be a level of the groups (`group.by`): ",
             paste(levels(tab$group), collapse = ", "), call. = FALSE)
    }
    has_facet <- "facet" %in% names(tab)
    seg <- data.frame(bar = tab$condition, group = tab$group)
    if (has_facet) seg <- cbind(facet = tab$facet, seg)
    seg$count <- tab$count
    seg$total <- tab$total
    seg$prop <- tab$prop
    seg$value <- if (value == "count") tab$count else tab$prop
    panel <- if (has_facet) droplevels(seg$facet) else factor(rep("all", nrow(seg)))

    # Positions: the bars of each panel, after those of the panels before
    seg$pos <- NA_real_
    offset <- 0
    for (pl in levels(panel)) {
        rows <- panel == pl
        bars <- levels(droplevels(seg$bar[rows]))
        if (!is.null(order.by)) {
            in_group <- rows & seg$group == order.by
            share <- seg$prop[in_group][match(bars, as.character(seg$bar[in_group]))]
            # order() is stable: ties keep the order of the levels
            bars <- bars[order(-share)]
        }
        seg$pos[rows] <- offset + match(as.character(seg$bar[rows]), bars)
        offset <- offset + length(bars)
    }

    # Stack each bar in the order of the groups: lower end of a segment = sum
    # of the segments stacked before it
    seg$lower <- NA_real_
    for (idx in split(seq_len(nrow(seg)), seg$pos)) {
        idx <- idx[order(seg$group[idx])]
        v <- seg$value[idx]
        before <- if (first.at.start) cumsum(v) - v else rev(cumsum(rev(v))) - v
        seg$lower[idx] <- before
    }
    seg$upper <- seg$lower + seg$value
    seg <- seg[order(seg$pos, seg$group), , drop = FALSE]
    rownames(seg) <- NULL
    cols <- c(
        if (has_facet) "facet",
        "bar",
        "pos",
        "group",
        "count",
        "total",
        "prop",
        "value",
        "lower",
        "upper"
    )
    seg[cols]
}


#' Bands joining the segments of neighbouring bars
#'
#' For each pair of neighbouring bars of a panel and each group, a polygon
#' from the group's segment on the first bar (its edge facing the second) to
#' its segment on the second, in a smooth S ("curve", as in alluvial plots)
#' or a straight line. Groups empty in both bars get no band; bars of
#' different panels are never joined.
#'
#' @param seg Result of \code{composition_segments()}
#'
#' @param half Half the width of a bar, in units of the category axis
#'
#' @param shape "curve" or "straight"
#'
#' @param n.points Points along each curved edge
#'
#' @return A data frame: `facet` (if any), `group`, `id` (one per polygon),
#'         `pos` and `val` (the polygon's corners). No rows if there is
#'         nothing to join
#'
#' @keywords internal
composition_bands <- function(seg,
                              half,
                              shape = c("curve", "straight"),
                              n.points = 40) {
    shape <- match.arg(shape)
    has_facet <- "facet" %in% names(seg)
    n <- if (shape == "curve") n.points else 2
    t <- seq(0, 1, length.out = n)
    # Smoothstep keeps the band level where it meets each bar
    s <- if (shape == "curve") 3 * t^2 - 2 * t^3 else t
    panel <- if (has_facet) droplevels(seg$facet) else factor(rep("all", nrow(seg)))
    bands <- list()
    id <- 0L
    for (pl in levels(panel)) {
        p_seg <- seg[panel == pl, , drop = FALSE]
        pos <- sort(unique(p_seg$pos))
        for (k in seq_len(length(pos) - 1)) {
            a <- p_seg[p_seg$pos == pos[k], , drop = FALSE]
            b <- p_seg[p_seg$pos == pos[k + 1], , drop = FALSE]
            b <- b[match(a$group, b$group), , drop = FALSE]
            x <- (pos[k] + half) + (pos[k + 1] - pos[k] - 2 * half) * t
            for (g in seq_len(nrow(a))) {
                if (a$value[g] == 0 && b$value[g] == 0) next
                top <- a$upper[g] + (b$upper[g] - a$upper[g]) * s
                bottom <- a$lower[g] + (b$lower[g] - a$lower[g]) * s
                id <- id + 1L
                # Along the top edge, then back along the bottom one
                band <- data.frame(
                    group = a$group[g],
                    id = id,
                    pos = c(x, rev(x)),
                    val = c(top, rev(bottom))
                )
                if (has_facet) band <- cbind(facet = a$facet[g], band)
                bands[[id]] <- band
            }
        }
    }
    if (length(bands) == 0) {
        empty <- data.frame(
            group = seg$group[0],
            id = integer(0),
            pos = numeric(0),
            val = numeric(0)
        )
        if (has_facet) empty <- cbind(facet = seg$facet[0], empty)
        return(empty)
    }
    do.call(rbind, bands)
}
