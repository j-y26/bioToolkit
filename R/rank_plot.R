# Rank plot of any rankable metric (addition 2026-11-16)
#
# Design
# ------
# `plot_rank()` draws features (genes, regulons, metabolites, gene sets...)
# along their rank on x, with a value on y: usually the metric they are
# ranked by, e.g. the regulon specificity score (RSS) of SCENIC, the reporter
# z-score of reporter metabolites, a log2 fold change or an NES. It
# generalizes the user's reporter-metabolite rank plot: any table, any two
# numeric columns (they can be the same).
#
# * Ties (user decision 2026-11-16, "option c"). Every point gets its own x
#   position, so tied points never sit on top of each other. Ties of
#   `rank.by` are broken by `value`, then by the feature name (radix sort, so
#   independent of the locale), never by the row order of the input: the
#   same data in another order gives the same plot and labels the same
#   features. The source data hold both the position (`rank`) and the tied,
#   competition rank (`rank_tied`, "1, 2, 2, 4") to report. `label.n` counts
#   positions; when the cutoff falls inside a run of ties, a message names
#   the tied features left out (labelling them all would flood the plot when
#   a metric has many ties, e.g. many zeros).
# * Which end is labelled (`label.side = "auto"`, user decision): both ends
#   for a signed metric (values on both sides of 0: z-scores, fold changes,
#   NES), where both directions mean something; only the first ranks
#   otherwise (RSS, specificity, -log10 p), where the end of the ranking is
#   not of interest. "first" / "last" name ends of the ranking rather than
#   high / low values, so that they keep their meaning with
#   `decreasing = FALSE` (e.g. p-values ranked from the smallest).
# * Labels (user decision 2026-11-16): each label is placed on its own
#   next to its point (no layout in columns), scattered around the corners
#   the curve leaves empty. ggrepel was tried first and dropped: in small
#   panels with long names (e.g. Hallmark sets) it left labels on top of
#   each other, piled against the panel edges. The labels are placed here
#   instead (`rank_place_labels()`), in inches, from the panel size: each
#   label goes to the nearest position around its point (centres of a
#   0.025 in grid over the panel) that is inside the panel and clear of the other
#   labels, of every point and of the leader lines (no line crosses a label
#   or another line). Positions on the side of the curve cost more, so that
#   labels of the first ranks sit towards higher ranks and those of the
#   last ranks towards lower ranks, where a monotone curve leaves room.
#   Labels nearest the ends are placed first. When a panel is too crowded,
#   the conditions are relaxed (lines, then points, then labels) with a
#   warning. The layout is deterministic (no seed) and needs no ggrepel.
#   For it, the panel limits are set by the function (ranks +-4%, values
#   +-10%, drawn as invisible corners with unexpanded scales), so the panel
#   size must be set with `panel.width` / `panel.height`. Labels are small
#   (7 pt, under the 9 pt tick labels) with a thin leader line each.
# * Colours: by default light grey points, the labelled ones in red (the
#   volcano's "up" red) and 1.3 times larger. `color.by` colours by a
#   column: numeric values with reversed magma (user decision; trimmed to
#   0.05-0.92 of the ramp as in the user's own plot, so that the lowest
#   values are not white and the highest not black), discrete ones with
#   group colours. Points with larger values are drawn on top.
#   `color.limits` fixes the colour scale, so that several figures (e.g.
#   time points) share it; values beyond it get the end colours.
# * Size: small panels (2.5 x 1.8 in) for a small space in a figure; text
#   sizes are those of `theme_sizes()`, as everywhere. `split.by` ranks
#   within each level and draws one panel per level (free scales), e.g. the
#   RSS of each cell type.


#' @title Rank plot of a metric
#'
#' @description Draws features (genes, regulons, metabolites, gene sets...)
#'              along their rank on x, with a value on y: usually the metric
#'              they are ranked by, such as the regulon specificity score
#'              (RSS) of SCENIC, reporter metabolite z-scores, log2 fold
#'              changes or NES. The features at the ends of the ranking, or
#'              chosen ones, are named. Panels have a fixed size; save the
#'              plot with \code{\link{save_plot}}.
#'
#' @details
#' \strong{Ties.} Each feature gets its own position on x, so that tied
#' features are not drawn on top of each other. Ties of `rank.by` are broken
#' by `value`, then by the feature name (alphabetically), never by the row
#' order of `data`, so the same data in any order give the same plot. The
#' source data hold the position (`rank`) and the tied rank to report
#' (`rank_tied`: competition ranking, e.g. 1, 2, 2, 4). `label.n` counts
#' positions; when the cutoff falls inside a run of ties, a message names
#' the tied features left unlabelled.
#'
#' \strong{Labels.} Each label is placed next to its point, at the nearest
#' position clear of the other labels, of all points and of the leader
#' lines, preferring the side the curve leaves empty (towards higher ranks
#' for the first ranks, towards lower ranks for the last ones); a leader
#' line joins it to its point. The layout is computed for the panel size,
#' so set the size with `panel.width` / `panel.height` (not by adding
#' \code{\link{panel_size_theme}} afterwards). When a panel is too crowded,
#' the labels are made smaller (with the default `label.size`, down to 5
#' pt); if they still do not fit, a warning suggests a smaller `label.n` or a
#' larger panel.
#'
#' @param data A data frame with one row per feature (per feature and
#'             `split.by` level)
#'
#' @param rank.by Numeric column the features are ranked by
#'
#' @param value Numeric column drawn on y, default is `rank.by`
#'
#' @param feature.col Column of the feature names. If NULL (default), the
#'                    row names
#'
#' @param decreasing Whether rank 1 is the largest value of `rank.by`
#'                   (default, TRUE) or the smallest (FALSE, e.g. for
#'                   p-values)
#'
#' @param split.by Column whose levels are ranked separately, one panel
#'                 each (e.g. the cell type of regulon specificity scores).
#'                 Default is NULL (one ranking)
#'
#' @param ncol Number of panel columns with `split.by`. If NULL (default), up
#'             to 4
#'
#' @param color.by Column the points are coloured by: numeric (a colour bar,
#'                 e.g. -log10 adjusted p-value) or discrete (groups). If
#'                 NULL (default), the points are grey and the labelled ones
#'                 red
#'
#' @param cols Colours. With a numeric `color.by`: a sequential palette name
#'             or colours from low to high; if NULL, reversed magma. With a
#'             discrete `color.by`: a palette name, or colours named by
#'             level or in level order. Without `color.by`: two colours,
#'             the labelled and the other points (or named `labelled`,
#'             `other`)
#'
#' @param color.limits Limits of a numeric colour scale, e.g. `c(0, 18)`, so
#'                     that several plots share one scale; values beyond
#'                     them get the end colours. If NULL (default), the range
#'                     of the data
#'
#' @param label Which features to name: "top" (default; the `label.n` first
#'              and / or last ranks, see `label.side`, plus any
#'              `label.features`), "selected" (only `label.features`) or
#'              "none"
#'
#' @param label.n Number of features named at each labelled end of the
#'                ranking (in each panel), default is 5
#'
#' @param label.side Ends of the ranking named with `label = "top"`: "auto"
#'                   (default; both ends when `rank.by` has values on both
#'                   sides of 0, such as z-scores or fold changes, the first
#'                   ranks otherwise), "first", "last" or "both"
#'
#' @param label.features Features to name, matched exactly
#'
#' @param label.size Size of the labels in points. If NULL (default), 7 pt,
#'                   or smaller (down to 5 pt, the same in every panel) when
#'                   the labels do not fit clear of each other at 7 pt
#'
#' @param label.wrap Labels longer than this many characters are wrapped
#'                   after spaces, hyphens, underscores or commas, default is
#'                   25; NULL or Inf for no wrapping
#'
#' @param italic.genes Whether the labels are set in italics, default is
#'                     FALSE. Set it to TRUE when the features are genes
#'
#' @param pt.size Size of the points, default is 1; the labelled ones are
#'                1.3 times larger
#'
#' @param alpha Opacity of the points, default is 1
#'
#' @param raster Whether to rasterize the points (needs ggrastr). If NULL
#'               (default), above 5,000 points
#'
#' @param raster.dpi Resolution of rasterized points, default is 1000
#'
#' @param title Plot title, default is NULL
#'
#' @param x.title Title of the x-axis, default is "Rank"
#'
#' @param y.title Title of the y-axis. If NULL (default), the name of
#'                `value`
#'
#' @param legend.title Title of the colour legend. If NULL (default), the
#'                     name of `color.by`
#'
#' @param panel.width Width of a panel in inches, default is 2.5
#'
#' @param panel.height Height of a panel in inches, default is 1.8
#'
#' @param legend.position Position of the legend, default is "right"
#'
#' @param bg Background colour, default is "transparent"
#'
#' @return A ggplot object; its source data
#'         (\code{\link{export_source_data}}) hold one row per feature:
#'         `feature`, `split` (with `split.by`), `rank`, `rank_tied`,
#'         `rank_value`, `value`, `color` (with `color.by`) and `labelled`
#'
#' @examples
#' set.seed(1)
#' # Reporter-metabolite-like z-scores: a signed metric, both ends named
#' df <- data.frame(
#'     metabolite = paste0("Metabolite_", 1:200),
#'     z = rnorm(200, sd = 1.5),
#'     padj = runif(200)
#' )
#' df$neg_log10_padj <- -log10(df$padj)
#' plot_rank(
#'     df,
#'     rank.by = "z",
#'     feature.col = "metabolite",
#'     color.by = "neg_log10_padj",
#'     y.title = "Reporter z-score"
#' )
#'
#' # Specificity scores per cell type: one panel each, first ranks named
#' rss <- expand.grid(regulon = paste0("TF", 1:60), cell_type = c("B", "T"))
#' rss$rss <- runif(nrow(rss), 0.1, 0.6)
#' plot_rank(rss, "rss", feature.col = "regulon", split.by = "cell_type", y.title = "RSS")
#'
#' @seealso \code{\link{plot_gsea_rank}}, \code{\link{plot_volcano}},
#'          \code{\link{save_plot}}
#'
#' @import ggplot2
#' @export
#'
plot_rank <- function(
    data,
    rank.by,
    value = rank.by,
    feature.col = NULL,
    decreasing = TRUE,
    split.by = NULL,
    ncol = NULL,
    color.by = NULL,
    cols = NULL,
    color.limits = NULL,
    label = c("top", "selected", "none"),
    label.n = 5,
    label.side = c("auto", "first", "last", "both"),
    label.features = NULL,
    label.size = NULL,
    label.wrap = 25,
    italic.genes = FALSE,
    pt.size = 1,
    alpha = 1,
    raster = NULL,
    raster.dpi = 1000,
    title = NULL,
    x.title = "Rank",
    y.title = NULL,
    legend.title = NULL,
    panel.width = 2.5,
    panel.height = 1.8,
    legend.position = "right",
    bg = "transparent"
) {
    label <- match.arg(label)
    label.side <- match.arg(label.side)
    if (!is.data.frame(data)) stop("`data` must be a data frame", call. = FALSE)
    if (!(is.numeric(label.n) && length(label.n) == 1 && !is.na(label.n) && label.n >= 0)) {
        stop("`label.n` must be a single number >= 0", call. = FALSE)
    }
    if (label == "selected" && length(label.features) == 0) {
        stop("`label = \"selected\"` needs `label.features`", call. = FALSE)
    }

    # Tidy table: one row per feature, with the ranking and plotted columns
    df <- rank_input(
        data,
        rank.by = rank.by,
        value = value,
        feature.col = feature.col,
        split.by = split.by,
        color.by = color.by
    )

    # Positions and tied ranks, within each split.by level
    df <- rank_positions(df, decreasing)

    # Ends named: both for a signed metric, the first ranks otherwise
    if (label.side == "auto") {
        signed <- any(df$rank_value < 0) && any(df$rank_value > 0)
        label.side <- if (signed) "both" else "first"
    }
    df$labelled <- rank_label_flags(df, label, label.n, label.side, label.features)

    sz <- theme_sizes()
    lw <- pt_to_linewidth(sz$line)

    # Colours: by a column, or labelled vs other points
    numeric_colour <- !is.null(color.by) && is.numeric(df$color)
    if (is.null(color.by)) {
        point_cols <- rank_highlight_colours(cols)
        df$.fill <- ifelse(df$labelled, "labelled", "other")
        # Labelled points on top
        df <- df[order(df$labelled), , drop = FALSE]
        fill_scale <- scale_fill_manual(
            values = point_cols,
            guide = "none"
        )
    } else if (numeric_colour) {
        pal <- if (is.null(cols)) {
            # Reversed magma, trimmed as in the user's reporter-metabolite
            # plot: neither white nor black at the ends
            viridisLite::magma(256, begin = 0.05, end = 0.92, direction = -1)
        } else {
            resolve_continuous_palette(cols, diverging = FALSE)
        }
        df$.fill <- df$color
        # Larger values on top; labelled points last of all
        df <- df[order(df$labelled, df$color, na.last = FALSE), , drop = FALSE]
        fill_scale <- scale_fill_gradientn(
            colours = pal,
            limits = color.limits,
            breaks = scales::breaks_extended(4),
            oob = scales::squish,
            na.value = "grey85",
            name = bold_title(if (is.null(legend.title)) color.by else legend.title),
            guide = guide_colourbar(
                theme = theme(
                    legend.key.height = unit(min(0.9, 0.6 * panel.height), "in"),
                    legend.key.width = unit(0.12, "in")
                )
            )
        )
    } else {
        df$color <- as.factor(df$color)
        df$.fill <- df$color
        df <- df[order(df$labelled), , drop = FALSE]
        fill_scale <- scale_fill_manual(
            values = group_colours(cols, df$color)$cols,
            name = bold_title(if (is.null(legend.title)) color.by else legend.title),
            guide = guide_legend(override.aes = list(size = 2.5))
        )
    }
    rownames(df) <- NULL
    df$.size <- ifelse(df$labelled, pt.size * 1.3, pt.size)

    # Limits of each panel, set here (not by the scales) so that the labels
    # can be placed in inches: ranks 1 to n plus 4% on each side, values
    # plus 10% above and below (room for labels at the ends of the curve)
    lims <- rank_panel_limits(df)

    # Label positions, panel by panel
    labs_df <- NULL
    if (any(df$labelled)) {
        df$.label <- ifelse(df$labelled, wrap_name(df$feature, label.wrap), "")
        # Default size: 7 pt, or smaller (to 5 pt) when the labels do not
        # fit at 7 pt; a size given is kept
        sizes <- if (is.null(label.size)) seq(7, 5, by = -0.5) else label.size
        labs_df <- rank_label_positions(
            df,
            lims = lims,
            panel.width = panel.width,
            panel.height = panel.height,
            sizes = sizes,
            italic = italic.genes,
            decreasing = decreasing
        )
        label.size <- attr(labs_df, "label.size")
    }

    layer_fun <- points_layer_fun(raster, nrow(df), 5000, raster.dpi)
    p <- ggplot(df, aes(x = .data$rank, y = .data$value)) +
        # Invisible corners of each panel: the scales span exactly the
        # limits the labels were placed in
        geom_blank(
            data = lims,
            aes(x = .data$x, y = .data$y)
        )
    # Zero line when the values cross 0 (in any panel)
    if (any(df$value < 0) && any(df$value > 0)) {
        p <- p + geom_hline(
            yintercept = 0,
            colour = "black",
            linewidth = lw
        )
    }
    # Leader lines under the points, so that each starts at its point's edge
    if (!is.null(labs_df)) {
        p <- p + geom_segment(
            data = labs_df,
            aes(
                x = .data$rank,
                y = .data$value,
                xend = .data$seg_x,
                yend = .data$seg_y
            ),
            colour = "grey30",
            linewidth = lw
        )
    }
    p <- p +
        layer_fun(geom_point(
            aes(fill = .data$.fill, size = .data$.size),
            shape = 21,
            colour = "black",
            stroke = pt_to_stroke(0.15),
            alpha = alpha
        )) +
        scale_size_identity() +
        fill_scale +
        scale_x_continuous(
            breaks = scales::breaks_extended(4),
            expand = expansion(0)
        ) +
        scale_y_continuous(expand = expansion(0)) +
        labs(
            title = title,
            x = x.title,
            y = if (is.null(y.title)) value else y.title
        ) +
        theme_bioToolkit("axis", bg = bg) +
        theme(legend.position = legend.position)
    if (!is.null(labs_df)) {
        p <- p + geom_text(
            data = labs_df,
            aes(
                x = .data$label_x,
                y = .data$label_y,
                label = .data$.label
            ),
            size = pt_to_size(label.size),
            fontface = if (italic.genes) "italic" else "plain",
            lineheight = RANK_LABEL_LINEHEIGHT,
            colour = "black"
        )
    }

    # One panel per split.by level, each with its own ranks and values
    if (!is.null(split.by)) {
        n_panels <- nlevels(df$split)
        layout <- panel_layout(n_panels, ncol, by.row = TRUE, default_ncol = min(n_panels, 4L))
        p <- p +
            facet_wrap(
                vars(.data$split),
                ncol = layout$ncol,
                scales = "free"
            ) +
            panel_size_theme(
                panel.width,
                panel.height,
                ncol = layout$ncol,
                nrow = layout$nrow
            )
    } else {
        p <- p + panel_size_theme(panel.width, panel.height)
    }

    src <- df[setdiff(names(df), c(".fill", ".size", ".label"))]
    if (is.null(split.by)) src$split <- NULL
    if (is.null(color.by)) src$color <- NULL
    attr(p, "source_data") <- src
    p
}


#' Tidy table of a rank plot
#'
#' @param data Input data frame
#'
#' @param rank.by,value,feature.col,split.by,color.by Columns, as in
#'        \code{\link{plot_rank}}
#'
#' @return A data frame with `feature`, `split` (factor, one level without
#'         `split.by`), `rank_value`, `value` and `color` (NULL without
#'         `color.by`), without the rows missing `rank.by` or `value`
#'
#' @keywords internal
rank_input <- function(data,
                       rank.by,
                       value,
                       feature.col,
                       split.by,
                       color.by) {
    check_col <- function(col, arg, numeric = FALSE) {
        if (is.null(col)) return(invisible())
        if (!(is.character(col) && length(col) == 1 && col %in% names(data))) {
            stop("`", arg, "` must be the name of a column of `data`", call. = FALSE)
        }
        if (numeric && !is.numeric(data[[col]])) {
            stop("`", arg, "` ('", col, "') must be a numeric column", call. = FALSE)
        }
    }
    check_col(rank.by, "rank.by", numeric = TRUE)
    check_col(value, "value", numeric = TRUE)
    check_col(feature.col, "feature.col")
    check_col(split.by, "split.by")
    check_col(color.by, "color.by")

    # Feature names: a column, or the row names (not automatic ones, which
    # would name features 1, 2, 3...)
    if (is.null(feature.col)) {
        if (.row_names_info(data) <= 0) {
            stop("`data` has no row names; give the column of feature names in `feature.col`", call. = FALSE)
        }
        feature <- rownames(data)
    } else {
        feature <- as.character(data[[feature.col]])
    }
    df <- data.frame(
        feature = feature,
        split = if (is.null(split.by)) factor("all") else droplevels(as.factor(data[[split.by]])),
        rank_value = data[[rank.by]],
        value = data[[value]],
        stringsAsFactors = FALSE
    )
    if (!is.null(color.by)) df$color <- data[[color.by]]

    # Rows that cannot be placed (NA / NaN / infinite rank or value) are left
    # out, with a count
    keep <- is.finite(df$rank_value) & is.finite(df$value)
    if (!is.null(split.by)) keep <- keep & !is.na(df$split)
    if (any(!keep)) {
        message(sum(!keep), " row(s) with a missing or infinite `rank.by` / `value`",
                if (!is.null(split.by)) " / `split.by`", " not drawn")
        df <- df[keep, , drop = FALSE]
    }
    if (nrow(df) == 0) stop("No rows left to rank", call. = FALSE)
    df$split <- droplevels(df$split)
    df
}


#' Rank positions and tied ranks, within each split
#'
#' Ties of `rank_value` are broken by `value`, then by the feature name
#' (radix sort, independent of the locale), so that the result does not
#' depend on the row order.
#'
#' @param df Table from \code{rank_input()}
#'
#' @param decreasing Whether rank 1 is the largest value
#'
#' @return `df` ordered by split and rank, with `rank` (position, 1 to n)
#'         and `rank_tied` (competition rank: tied values share the smallest
#'         rank of the run)
#'
#' @keywords internal
rank_positions <- function(df,
                           decreasing) {
    sgn <- if (decreasing) -1 else 1
    ord <- order(
        df$split,
        sgn * df$rank_value,
        sgn * df$value,
        df$feature,
        method = "radix"
    )
    df <- df[ord, , drop = FALSE]
    df$rank <- stats::ave(seq_len(nrow(df)), df$split, FUN = seq_along)
    df$rank_tied <- stats::ave(sgn * df$rank_value, df$split, FUN = function(x) rank(x, ties.method = "min"))
    df$rank_tied <- as.integer(df$rank_tied)
    rownames(df) <- NULL
    df
}


#' Which features of a rank plot are named
#'
#' @param df Table from \code{rank_positions()}
#'
#' @param label "top", "selected" or "none"
#'
#' @param label.n Number of features per labelled end, in each split
#'
#' @param label.side "first", "last" or "both"
#'
#' @param label.features Features to name, or NULL
#'
#' @return A logical vector along the rows of `df`
#'
#' @keywords internal
rank_label_flags <- function(df,
                             label,
                             label.n,
                             label.side,
                             label.features) {
    flags <- rep(FALSE, nrow(df))
    if (label == "none") return(flags)
    if (label == "top" && label.n > 0) {
        n <- stats::ave(df$rank, df$split, FUN = length)
        first <- df$rank <= label.n
        last <- df$rank > n - label.n
        if (label.side %in% c("first", "both")) flags <- flags | first
        if (label.side %in% c("last", "both")) flags <- flags | last
        # Ties at the cutoff: tied features on the other side of it share a
        # tied rank with a labelled one but are not named; say which
        left_out <- !flags & stats::ave(flags, df$split, df$rank_tied, FUN = any) > 0
        if (any(left_out)) {
            shown <- utils::head(df$feature[left_out], 5)
            message(
                sum(left_out), " feature(s) tied with a labelled one at the `label.n` cutoff are not ",
                "labelled (ties are broken by `value`, then by name): ",
                paste(shown, collapse = ", "),
                if (sum(left_out) > 5) ", ..."
            )
        }
    }
    if (length(label.features)) {
        label.features <- as.character(label.features)
        missing <- setdiff(label.features, df$feature)
        if (length(missing)) {
            warning("`label.features` not found among the ranked features: ", paste(missing, collapse = ", "),
                    call. = FALSE)
        }
        flags <- flags | df$feature %in% label.features
    }
    flags
}


#' Colours of the labelled and other points of a rank plot
#'
#' @param cols NULL, two colours (labelled, other) or colours named
#'             `labelled` / `other`
#'
#' @return Two colours named `labelled` and `other`
#'
#' @keywords internal
rank_highlight_colours <- function(cols) {
    # The volcano's "up" red, and a light grey that leaves the curve visible
    out <- c(labelled = "#C91F28", other = "grey80")
    if (is.null(cols)) return(out)
    if (!all(is_colour(cols))) stop("`cols` must be colours", call. = FALSE)
    if (is.null(names(cols))) {
        if (length(cols) > 2) stop("Without `color.by`, `cols` is two colours: labelled, other", call. = FALSE)
        names(cols) <- names(out)[seq_along(cols)]
    }
    unknown <- setdiff(names(cols), names(out))
    if (length(unknown)) {
        stop("`cols` must be named 'labelled' and / or 'other', not: ", paste(unknown, collapse = ", "),
             call. = FALSE)
    }
    out[names(cols)] <- cols
    out
}


# Line height of multi-line labels
RANK_LABEL_LINEHEIGHT <- 0.85


#' Limits of each panel of a rank plot
#'
#' @param df Plotted table, with `split`, `rank` and `value`
#'
#' @return A data frame with two rows per split (the lower-left and
#'         upper-right corners): `split`, `x`, `y`
#'
#' @keywords internal
rank_panel_limits <- function(df) {
    out <- lapply(levels(df$split), function(s) {
        d <- df[df$split == s, , drop = FALSE]
        n <- max(d$rank)
        # Ranks 1..n plus 4% on each side (one rank on each side for n = 1)
        pad_x <- if (n > 1) 0.04 * (n - 1) else 1
        y <- range(d$value)
        # Values plus 10% above and below; 1 on each side of a constant
        pad_y <- if (diff(y) > 0) 0.1 * diff(y) else max(1, abs(y[1]) * 0.1)
        data.frame(
            split = factor(s, levels = levels(df$split)),
            x = c(1 - pad_x, n + pad_x),
            y = c(y[1] - pad_y, y[2] + pad_y)
        )
    })
    do.call(rbind, out)
}


#' Sizes of labels in inches
#'
#' Measured on a PDF device, plus 10% of the width (raster devices may use a
#' slightly wider font) and a small padding on each side.
#'
#' @param labels Labels (may contain line breaks)
#'
#' @param size Font size in points
#'
#' @param italic Whether the labels are italic
#'
#' @return A list of `w` and `h`, in inches
#'
#' @keywords internal
rank_label_extent <- function(labels,
                              size,
                              italic) {
    # Measure on a throwaway device, then return to the current one
    old_dev <- grDevices::dev.cur()
    grDevices::pdf(NULL)
    on.exit({
        grDevices::dev.off()
        if (old_dev > 1) grDevices::dev.set(old_dev)
    }, add = TRUE)
    gp <- grid::gpar(
        fontsize = size,
        fontface = if (italic) 3 else 1,
        lineheight = RANK_LABEL_LINEHEIGHT
    )
    w <- vapply(labels, function(s) {
        g <- grid::textGrob(s, gp = gp)
        grid::convertWidth(grid::grobWidth(g), "in", valueOnly = TRUE)
    }, numeric(1), USE.NAMES = FALSE)
    n_lines <- lengths(strsplit(labels, "\n", fixed = TRUE))
    h <- (1 + (n_lines - 1) * RANK_LABEL_LINEHEIGHT) * size / 72
    pad <- 0.012
    list(w = 1.1 * w + 2 * pad, h = h + 2 * pad)
}


#' Positions of the labels of a rank plot
#'
#' Each panel is laid out on its own, in inches, from its limits and the
#' panel size (see \code{rank_place_labels()}). The layout is tried at each
#' size of `sizes` in turn, until the labels of every panel fit clear of
#' the points and of each other, so that all panels share one label size.
#'
#' @param df Plotted table, with `split`, `rank`, `value`, `labelled`,
#'           `.size` and `.label`
#'
#' @param lims Panel limits from \code{rank_panel_limits()}
#'
#' @param panel.width,panel.height Panel size in inches
#'
#' @param sizes Font sizes of the labels to try, in points, largest first
#'
#' @param italic Whether the labels are italic
#'
#' @param decreasing Whether rank 1 is the largest value (the curve falls)
#'
#' @return The labelled rows of `df`, with `label_x`, `label_y` (centre of
#'         the label) and `seg_x`, `seg_y` (end of the leader line), in data
#'         units; attribute `label.size`, the size used
#'
#' @keywords internal
rank_label_positions <- function(df,
                                 lims,
                                 panel.width,
                                 panel.height,
                                 sizes,
                                 italic,
                                 decreasing) {
    # Panels in inches, computed once for all sizes
    panels <- list()
    for (s in levels(df$split)) {
        d <- df[df$split == s, , drop = FALSE]
        if (!any(d$labelled)) next
        lim <- lims[lims$split == s, , drop = FALSE]
        sx <- panel.width / diff(lim$x)
        sy <- panel.height / diff(lim$y)
        px <- (d$rank - lim$x[1]) * sx
        py <- (d$value - lim$y[1]) * sy
        # Radius of a point of ggplot2 size `s`: s * .pt points across
        pr <- d$.size * ggplot2::.pt / 72 / 2
        # Obstacles: one point per 0.015 in cell (whole gene lists would
        # otherwise make tens of thousands), and every labelled point
        cell <- paste(round(px / 0.015), round(py / 0.015))
        keep <- which(d$labelled | !duplicated(cell))
        lab <- which(d$labelled)
        n <- max(d$rank)
        # Each end of the ranking is a group: the first half (labels
        # preferably towards higher ranks, i.e. right) and the second half
        # (towards lower ranks, left)
        first <- d$rank[lab] <= n / 2
        side <- ifelse(first, 1, -1)
        # Vertically, away from the curve: above it at the first ranks of a
        # falling curve, below it at the last ranks (the other way round for
        # a rising curve)
        vside <- if (decreasing) side else -side
        panels[[s]] <- list(
            d = d[lab, , drop = FALSE],
            lim = lim,
            sx = sx,
            sy = sy,
            x = px[lab],
            y = py[lab],
            side = side,
            vside = vside,
            # Orders in which the labels are placed (see below)
            # Each end from its extreme inwards, or the other way round
            # (the four combinations; see below)
            orders = local({
                e <- pmin(d$rank[lab], n + 1 - d$rank[lab])
                unique(list(
                    order(e),
                    order(ifelse(first, e, -e)),
                    order(ifelse(first, -e, e)),
                    order(-e)
                ))
            }),
            obs_x = px[keep],
            obs_y = py[keep],
            obs_r = pr[keep],
            # Each label's own point among the obstacles
            self = match(lab, keep)
        )
    }

    for (size in sizes) {
        out <- list()
        crowded <- character(0)
        for (s in names(panels)) {
            pn <- panels[[s]]
            ext <- rank_label_extent(pn$d$.label, size, italic)
            # Labels are placed one after the other, so the first ones take
            # the best spots and can wall in the others (around a dense end
            # of the curve). The layout is made with each end placed from
            # its extreme inwards or outwards (four combinations) and the
            # best kept: the lowest total cost (with the penalties of relaxed
            # conditions and of labels out of order), labels over points or
            # labels counting heavily
            place <- function(ord) {
                res <- rank_place_labels(
                    x = pn$x,
                    y = pn$y,
                    w = ext$w,
                    h = ext$h,
                    side = pn$side,
                    vside = pn$vside,
                    order = ord,
                    obs_x = pn$obs_x,
                    obs_y = pn$obs_y,
                    obs_r = pn$obs_r,
                    self = pn$self,
                    panel.width = panel.width,
                    panel.height = panel.height
                )
                # Labels stacked in the opposite order to their points (in
                # y, within an end of the ranking) are hard to follow: each
                # such pair adds 0.15 in
                inv <- rank_label_inversions(pn$y, res$cy, pn$side)
                list(res = res, score = sum(res$cost) + 0.15 * inv + 10 * sum(res$level >= 4))
            }
            best <- NULL
            for (k in seq_along(pn$orders)) {
                cand <- place(pn$orders[[k]])
                # The first order (from the ends inwards) is kept unless
                # another is clearly (15%) better
                if (k > 1) cand$score <- cand$score * 1.15
                if (is.null(best) || cand$score < best$score) best <- cand
            }
            res <- best$res
            if (any(res$level >= 4)) crowded <- c(crowded, s)
            d <- pn$d
            d$label_x <- pn$lim$x[1] + res$cx / pn$sx
            d$label_y <- pn$lim$y[1] + res$cy / pn$sy
            d$seg_x <- pn$lim$x[1] + res$ex / pn$sx
            d$seg_y <- pn$lim$y[1] + res$ey / pn$sy
            out[[s]] <- d
        }
        if (length(crowded) == 0) break
    }
    if (length(sizes) > 1 && size < sizes[1]) {
        message("Labels set at ", size, " pt to fit the panel")
    }
    if (length(crowded)) {
        where <- if (identical(crowded, "all")) "" else paste0(" (", paste(crowded, collapse = ", "), ")")
        warning("Some labels could not be placed clear of the points or of each other", where,
                "; use a smaller `label.n` or `label.size`, or a larger panel", call. = FALSE)
    }
    res <- do.call(rbind, unname(out))
    attr(res, "label.size") <- size
    res
}


#' Place labels around their points without overlaps
#'
#' All in inches, origin at the lower left of the panel. Each label is a
#' box; its candidate positions are the centres of a 0.025 in grid over the
#' panel, and the cost of a position is the gap between the point and the
#' nearest edge of the box (at least 0.025 in). Labels are placed one at a time in `order`, each at the
#' cheapest candidate that is inside the panel, clear of the labels placed
#' and of all points, whose leader line passes through no label and no
#' other point and crosses no other line, and whose box covers no line.
#' The cost is the gap, 1.8 times higher towards the curve horizontally
#' (`side`) and again vertically (`vside`), so that labels sit close to
#' their points, in the corner the curve leaves empty; a box across the
#' horizontal strip beside a point still to be labelled (on its side, at its
#' height) costs
#' 0.1 in more per strip, so that labels leave each other a way out and
#' form a staircase rather than a wall. A line may cross
#' another line (level 2, cost + 0.15 in) or pass over points (level 3,
#' + 0.3 in) when that keeps it shorter: a short line over a point reads
#' better than a long one across the plot. A line never passes through a
#' label. When no candidate is clear of labels and points, labels may cover
#' points (4), and as a last resort overlap other labels (5: the candidate
#' overlapping them least).
#'
#' @param x,y Positions of the labelled points
#'
#' @param w,h Sizes of the labels
#'
#' @param side 1 to prefer the right of the point, -1 the left
#'
#' @param vside 1 to prefer above the point, -1 below
#'
#' @param order Order in which the labels are placed
#'
#' @param obs_x,obs_y,obs_r Points to keep clear of, and their radii
#'
#' @param self Index of each label's own point among the obstacles
#'
#' @param panel.width,panel.height Panel size
#'
#' @return A data frame along the labels: centre (`cx`, `cy`), end of the
#'         leader line on the box (`ex`, `ey`), the `level` of the
#'         conditions met (1 = all) and the `cost` of the position, with
#'         the penalty of its level
#'
#' @keywords internal
rank_place_labels <- function(x,
                              y,
                              w,
                              h,
                              side,
                              vside,
                              order,
                              obs_x,
                              obs_y,
                              obs_r,
                              self,
                              panel.width,
                              panel.height) {
    # Gaps: between labels, and between the panel edge and a label
    gap <- 0.02
    edge <- 0.01
    # Candidate centres: a 0.025 in grid over the panel
    step <- 0.025
    grid_c <- expand.grid(
        x = seq(step / 2, panel.width, by = step),
        y = seq(step / 2, panel.height, by = step)
    )
    # Points along a leader line checked against boxes and points
    t_seg <- seq(0.1, 1, length.out = 10)

    out <- data.frame(cx = x, cy = y, ex = x, ey = y, level = 1L, cost = 0)
    placed <- rep(FALSE, length(x))
    # Placed labels (left, right, bottom, top) and leader lines
    bl <- br <- bb <- bt <- numeric(0)
    sx1 <- sy1 <- sx2 <- sy2 <- numeric(0)

    # Relaxed conditions cost extra (in inches of line): a leader line may
    # cross another line (level 2) or pass over points (level 3) when that
    # keeps it this much shorter
    penalty <- c(0, 0.15, 0.3)

    for (i in order) {
        cx <- grid_c$x
        cy <- grid_c$y
        # End of the leader line: the point of the box nearest the point
        ex <- pmin(pmax(x[i], cx - w[i] / 2), cx + w[i] / 2)
        ey <- pmin(pmax(y[i], cy - h[i] / 2), cy + h[i] / 2)
        # Cost: the gap between the point and the box, higher towards the
        # curve (horizontally and vertically)
        r <- sqrt((ex - x[i])^2 + (ey - y[i])^2)
        ux <- (ex - x[i]) / pmax(r, 1e-9)
        uy <- (ey - y[i]) / pmax(r, 1e-9)
        cost <- r *
            ifelse(ux * side[i] < -0.2, 1.8, 1) *
            ifelse(uy * vside[i] < -0.2, 1.8, 1)
        # Leave room for the labels still to come: a box across the
        # horizontal strip beside a point not yet labelled (on that point's
        # side, at its height) would wall its label in; each such
        # strip covered adds 0.1 in, so labels step aside into a staircase
        for (j in which(!placed & side == side[i])) {
            if (j == i) next
            across <- cy - h[i] / 2 < y[j] + 0.005 & cy + h[i] / 2 > y[j] - 0.005
            beside <- if (side[j] > 0) cx + w[i] / 2 > x[j] + 0.025 else cx - w[i] / 2 < x[j] - 0.025
            cost <- cost + 0.1 * (across & beside)
        }
        # Inside the panel and at least 0.025 in from the point, cheapest
        # first
        ok <- cx - w[i] / 2 >= edge & cx + w[i] / 2 <= panel.width - edge &
            cy - h[i] / 2 >= edge & cy + h[i] / 2 <= panel.height - edge &
            r >= 0.025
        idx <- which(ok)
        idx <- idx[order(cost[idx])]
        # Obstacles other than this label's own point
        ox <- obs_x[-self[i]]
        oy <- obs_y[-self[i]]
        orad <- obs_r[-self[i]]

        # Area of overlap of the candidates' boxes with the placed labels
        overlap_area <- function(k) {
            a <- numeric(length(k))
            for (b in seq_along(bl)) {
                dx <- pmin(cx[k] + w[i] / 2 + gap, br[b]) - pmax(cx[k] - w[i] / 2 - gap, bl[b])
                dy <- pmin(cy[k] + h[i] / 2 + gap, bt[b]) - pmax(cy[k] - h[i] / 2 - gap, bb[b])
                a <- a + pmax(dx, 0) * pmax(dy, 0)
            }
            a
        }
        # Whether any obstacle lies within `hw` x `hh` (plus its radius) of
        # the positions (px, py)
        hits_points <- function(px, py, hw, hh) {
            if (length(ox) == 0) return(rep(FALSE, length(px)))
            dx <- abs(outer(px, ox, "-")) < outer(hw, orad, "+")
            dy <- abs(outer(py, oy, "-")) < outer(hh, orad, "+")
            rowSums(dx & dy) > 0
        }
        # Leader lines of candidates `k`, checked once: `clear` (through no
        # placed label, and no placed line under the box: required),
        # `no_pts` (over no other point), `no_cross` (crossing no placed
        # line)
        line_checks <- function(k) {
            nk <- length(k)
            nt <- length(t_seg)
            # Points along each line: one row per candidate, one column per
            # step along the line
            lx <- x[i] + outer(ex[k] - x[i], t_seg)
            ly <- y[i] + outer(ey[k] - y[i], t_seg)
            inside <- matrix(FALSE, nk, nt)
            for (b in seq_along(bl)) {
                inside <- inside | (lx > bl[b] & lx < br[b] & ly > bb[b] & ly < bt[b])
            }
            clear <- rowSums(inside) == 0
            # Other points only beyond 0.06 in of the own point: in a dense
            # stretch of the curve the neighbouring points touch it, and any
            # line leaving it would cross them
            far <- outer(sqrt((ex[k] - x[i])^2 + (ey[k] - y[i])^2), t_seg) > 0.06
            hit <- matrix(FALSE, nk, nt)
            if (any(far)) {
                hit[far] <- hits_points(lx[far], ly[far], rep(0, sum(far)), rep(0, sum(far)))
            }
            no_pts <- rowSums(hit) == 0
            no_cross <- rep(TRUE, nk)
            for (s in seq_along(sx1)) {
                # Placed line under the box
                px_s <- sx1[s] + t_seg * (sx2[s] - sx1[s])
                py_s <- sy1[s] + t_seg * (sy2[s] - sy1[s])
                under <- abs(outer(cx[k], px_s, "-")) < w[i] / 2 & abs(outer(cy[k], py_s, "-")) < h[i] / 2
                clear <- clear & rowSums(under) == 0
                # Crossings among the points (within 0.06 in of the own
                # point) do not show and are allowed
                no_cross <- no_cross & !segments_cross(x[i], y[i], ex[k], ey[k],
                                                       sx1[s], sy1[s], sx2[s], sy2[s],
                                                       min.dist = 0.06)
            }
            list(clear = clear, no_pts = no_pts, no_cross = no_cross)
        }
        # Level of each candidate: 1 (all conditions met), 2 (its line
        # crosses another), 3 (its line passes over points) or NA (its box
        # overlaps a label or a point, or its line passes through a label
        # or under its box lies a line)
        candidate_level <- function(k) {
            lev <- rep(NA_integer_, length(k))
            ok <- overlap_area(k) == 0
            ok[ok] <- !hits_points(cx[k[ok]], cy[k[ok]], rep(w[i] / 2, sum(ok)), rep(h[i] / 2, sum(ok)))
            if (!any(ok)) return(lev)
            lc <- line_checks(k[ok])
            lev[ok] <- ifelse(!lc$clear, NA_integer_,
                              ifelse(!lc$no_pts, 3L, ifelse(lc$no_cross, 1L, 2L)))
            lev
        }

        # Cheapest candidate, its cost raised by the penalty of its level;
        # candidates are scanned in order of cost, until none left can beat
        # the best found
        chosen <- NA_integer_
        level <- NA_integer_
        best <- Inf
        for (start in seq_len(ceiling(length(idx) / 256))) {
            k <- idx[((start - 1) * 256 + 1):min(length(idx), start * 256)]
            if (cost[k[1]] >= best) break
            lev <- candidate_level(k)
            score <- cost[k] + penalty[lev]
            if (any(!is.na(score)) && min(score, na.rm = TRUE) < best) {
                j <- which.min(score)
                best <- score[j]
                chosen <- k[j]
                level <- lev[j]
            }
        }
        if (is.na(chosen)) {
            # Labels may cover points (level 4), or as a last resort overlap
            # other labels (5: the candidate overlapping them least)
            free <- if (length(idx)) idx[overlap_area(idx) == 0] else integer(0)
            if (length(free)) {
                chosen <- free[1]
                level <- 4L
            } else if (length(idx)) {
                chosen <- idx[which.min(overlap_area(idx) + 1e-6 * cost[idx])]
                level <- 5L
            } else {
                # No candidate inside the panel (label wider than the panel)
                chosen <- which.min(cost)
                level <- 5L
            }
        }
        best_score <- if (is.finite(best)) best else cost[chosen]

        out$cx[i] <- cx[chosen]
        out$cy[i] <- cy[chosen]
        out$ex[i] <- ex[chosen]
        out$ey[i] <- ey[chosen]
        out$level[i] <- level
        out$cost[i] <- best_score
        placed[i] <- TRUE
        bl <- c(bl, cx[chosen] - w[i] / 2)
        br <- c(br, cx[chosen] + w[i] / 2)
        bb <- c(bb, cy[chosen] - h[i] / 2)
        bt <- c(bt, cy[chosen] + h[i] / 2)
        sx1 <- c(sx1, x[i])
        sy1 <- c(sy1, y[i])
        sx2 <- c(sx2, ex[chosen])
        sy2 <- c(sy2, ey[chosen])
    }
    out
}


#' Number of pairs of labels stacked in the opposite order to their points
#'
#' @param y Heights of the points
#'
#' @param cy Heights of the labels' centres
#'
#' @param group Group of each label (an end of the ranking); pairs are
#'              counted within groups
#'
#' @return The number of pairs whose labels and points are in opposite
#'         vertical order
#'
#' @keywords internal
rank_label_inversions <- function(y,
                                  cy,
                                  group) {
    n_inv <- 0
    for (g in unique(group)) {
        k <- which(group == g)
        if (length(k) < 2) next
        dy <- outer(y[k], y[k], "-")
        dc <- outer(cy[k], cy[k], "-")
        n_inv <- n_inv + sum((dy * dc < 0)[upper.tri(dy)])
    }
    n_inv
}


#' Whether segments cross
#'
#' Proper crossings only: segments that touch at an end (e.g. two leader
#' lines from neighbouring points) do not count, nor do crossings within
#' `min.dist` of the start of the first segment.
#'
#' @param ax1,ay1,ax2,ay2 First segment(s), vectorized
#'
#' @param bx1,by1,bx2,by2 Second segment
#'
#' @param min.dist Crossings nearer than this to (ax1, ay1) are ignored
#'
#' @return A logical vector along the first segments
#'
#' @keywords internal
segments_cross <- function(ax1,
                           ay1,
                           ax2,
                           ay2,
                           bx1,
                           by1,
                           bx2,
                           by2,
                           min.dist = 0) {
    orient <- function(px, py, qx, qy, rx, ry) (qx - px) * (ry - py) - (qy - py) * (rx - px)
    d1 <- orient(bx1, by1, bx2, by2, ax1, ay1)
    d2 <- orient(bx1, by1, bx2, by2, ax2, ay2)
    d3 <- orient(ax1, ay1, ax2, ay2, bx1, by1)
    d4 <- orient(ax1, ay1, ax2, ay2, bx2, by2)
    cross <- (d1 * d2 < 0) & (d3 * d4 < 0)
    # Distance of the crossing from the start of the first segment
    t <- d1 / (d1 - d2)
    dist <- t * sqrt((ax2 - ax1)^2 + (ay2 - ay1)^2)
    cross & !is.na(dist) & dist > min.dist
}


#' Wrap long names after spaces, hyphens, underscores or commas
#'
#' @param x Names
#'
#' @param width Maximum characters per line; NULL or Inf for no wrapping
#'
#' @return `x`, with line breaks; a hyphen stays at the end of its line
#'
#' @keywords internal
wrap_name <- function(x,
                      width) {
    if (is.null(width) || !is.finite(width)) return(x)
    vapply(x, function(s) {
        if (is.na(s) || nchar(s) <= width) return(s)
        # Pieces ending with a breakable character
        pieces <- regmatches(s, gregexpr("[^ _,-]*[ _,-]*", s))[[1]]
        pieces <- pieces[nzchar(pieces)]
        lines <- character(0)
        cur <- ""
        for (pc in pieces) {
            if (nchar(cur) > 0 && nchar(trimws(paste0(cur, pc))) > width) {
                lines <- c(lines, trimws(cur))
                cur <- pc
            } else {
                cur <- paste0(cur, pc)
            }
        }
        paste(c(lines, trimws(cur)), collapse = "\n")
    }, character(1), USE.NAMES = FALSE)
}

