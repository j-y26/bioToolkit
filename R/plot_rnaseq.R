# Volcano plot of differential expression results (roadmap item R2)
#
# Design
# ------
# * One function, `plot_volcano()`, replaces `plot_volcano()` and
#   `plot_volcano2()` of bioToolkit <= 0.5.0. The two only differed in how the
#   numbers of up- and down-regulated features were shown, which is now
#   `annotate.counts = "corner"` (old plot_volcano) or `"arrow"` (old
#   plot_volcano2).
# * DE results are read with `standardize_de()` (F5), so edgeR, DESeq2, limma
#   and Seurat tables work without naming columns, p-values of 0 sit just above
#   the others instead of at 1e-300, and NA adjusted p-values (DESeq2
#   independent filtering) are handled in one place.
# * The y-axis is the -log10 *adjusted* p-value (decided 2026-10-12). The
#   colours are decided by the same adjusted p-value, so the dashed line at
#   `padj.cutoff` separates exactly the coloured points from the grey ones and
#   the figure can be read without the legend.
# * Every point is drawn at its own value: no capping of extreme logFC or
#   p-values (decided 2026-10-12).
# * Labels are bold (bold italic for genes), to stay legible over the points
#   and at print sizes (decided 2026-10-13).
# * Labels are chosen explicitly (`label` = "top", "selected" or "none"). The
#   old functions silently replaced the requested genes by the top 10 per
#   direction when many were requested; that behaviour is gone.
# * Non-significant points are drawn first, so that they never hide DE points.
# * Fixed panel size (F2) and `theme_bioToolkit("axis")` (F3). With more than
#   5,000 features the points are rasterized (`raster.dpi`, default 1000),
#   while labels, lines and text stay vector.


#' @title Volcano plot of differential expression results
#'
#' @description Plots the log2 fold change of each feature against its
#'              -log10 adjusted p-value, colouring features that are
#'              significantly up- or down-regulated, with optional labels and
#'              counts of DE features. Reads edgeR, DESeq2, limma and Seurat
#'              results directly (see \code{\link{standardize_de}}).
#'
#' @param data DE results of one comparison: a data frame (or a DESeq2
#'             `DESeqResults` / edgeR `TopTags`), with feature names as row
#'             names or in `feature.col`
#'
#' @param method The tool that produced the results: "auto" (default) detects
#'               it from the column names; or one of "edgeR", "DESeq2",
#'               "limma" and "seurat"
#'
#' @param feature.col Column of feature (gene) names. If NULL (default), the
#'                    row names (or the `gene` column of Seurat results)
#'
#' @param logfc.col Column of the log2 fold change. If NULL (default), taken
#'                  from `method`. To plot a shrunken log fold change (DESeq2
#'                  \code{lfcShrink()}, e.g. with apeglm), which is less noisy
#'                  for lowly expressed genes, give its column here
#'
#' @param p.col Column of the raw p-value. If NULL (default), taken from
#'              `method`. Only used to compute adjusted p-values
#'              (Benjamini-Hochberg) when the table has none
#'
#' @param padj.col Column of the adjusted p-value. If NULL (default), taken
#'                 from `method`. It is both plotted on the y-axis and used to
#'                 decide significance
#'
#' @param padj.cutoff Adjusted p-value below which a feature is significant,
#'                    default is 0.05. Drawn as a dashed horizontal line
#'
#' @param logfc.cutoff Absolute log2 fold change a significant feature must
#'                     exceed, default is 0 (any change). When above 0, drawn
#'                     as dashed vertical lines
#'
#' @param label Which features to label:
#'              \describe{
#'                \item{"top"}{(default) the `label.n` most significant
#'                      features of each direction, ranked by `label.rank`,
#'                      plus any `label.features`}
#'                \item{"selected"}{only `label.features`, whether or not
#'                      they are significant}
#'                \item{"none"}{no labels}
#'              }
#'
#' @param label.n Number of features to label per direction with
#'                `label = "top"`, default is 10. Every requested label is
#'                drawn; when long labels (e.g. Ensembl IDs) crowd the same
#'                corner they can overlap, so lower `label.n` or widen
#'                `panel.width`
#'
#' @param label.features Features to label. Always labelled with
#'                       `label = "selected"`, and in addition to the top
#'                       features with `label = "top"`. Matched exactly, so
#'                       the case of gene symbols matters (e.g. "CD3E" in
#'                       human, "Cd3e" in mouse and rat)
#'
#' @param label.rank How the top features are ranked with `label = "top"`:
#'                   "padj" (default; smallest adjusted p-value, ties broken
#'                   by the largest absolute log2 fold change), "logfc"
#'                   (largest absolute log2 fold change) or "product"
#'                   (-log10 adjusted p-value x absolute log2 fold change,
#'                   which favours features that are both)
#'
#' @param label.size Size of the labels in mm. If NULL (default), the small
#'                   small text size of the theme (9 pt)
#'
#' @param italic.genes Whether labels are gene names, set in italics. Default
#'                     is TRUE; set FALSE for features that are not genes
#'                     (e.g. peaks or pathways). Labels are bold either way
#'
#' @param annotate.counts How to show the numbers of up- and down-regulated
#'                        features: "corner" (default; in the bottom corners),
#'                        "arrow" (arrows naming the direction at the top of
#'                        the panel, with the counts below them) or "none"
#'
#' @param symmetric Whether the x-axis is symmetric around 0, default is TRUE
#'
#' @param cols Colours of the up-regulated, down-regulated and not significant
#'             features: a vector named `up`, `down` and `ns` (any of them;
#'             the others keep their default), or three colours in that
#'             order. NULL (default) is
#'             `c(up = "firebrick4", down = "midnightblue", ns = "grey75")`
#'
#' @param up.label,down.label,ns.label Names of the three classes in the
#'        legend (and next to the arrows with `annotate.counts = "arrow"`).
#'        Default "Up-regulated", "Down-regulated" and "Not significant"
#'
#' @param pt.size Point size, default is 1
#'
#' @param alpha Point opacity, default is 0.6
#'
#' @param raster Whether to rasterize the points while keeping labels, lines
#'               and text as vectors: NULL (default) rasterizes when more than
#'               5,000 features are drawn, or TRUE / FALSE
#'
#' @param raster.dpi Resolution of the rasterized points, default is 1000.
#'                   Match the `dpi` used when saving raster (e.g. png) output
#'
#' @param panel.width Width of the panel in inches, default is 3
#'
#' @param panel.height Height of the panel in inches. Default is
#'                     `panel.width`, a square panel, so that setting
#'                     `panel.width` alone resizes the plot
#'
#' @param legend.position Position of the legend. If NULL (default), "right",
#'                        or "none" with `annotate.counts = "arrow"` (the
#'                        arrows already name the directions). A legend at the
#'                        "top" is usually wider than the panel, which widens
#'                        the saved figure
#'
#' @param bg Background colour, default is "transparent"
#'
#' @param seed Random seed for the placement of labels, default is 42
#'
#' @details
#' A feature is up-regulated when its adjusted p-value is below `padj.cutoff`
#' and its log2 fold change is above `logfc.cutoff`, down-regulated when below
#' `-logfc.cutoff`, and not significant otherwise. Features without a log fold
#' change or an adjusted p-value are not drawn.
#'
#' The plot has a fixed panel size; save it with \code{\link{save_plot}},
#' which sizes the file to fit. Its data (one row per drawn feature, with
#' columns `direction`, `neg_log10_padj` and `labelled` added to the
#' standardized DE table) can be written with
#' \code{\link{export_source_data}}.
#'
#' @return A ggplot object
#'
#' @examples
#' set.seed(1)
#' de <- data.frame(
#'     logFC = rnorm(2000, sd = 1.5),
#'     logCPM = runif(2000, 0, 10),
#'     PValue = runif(2000)^4,
#'     row.names = paste0("Gene", 1:2000)
#' )
#' de$FDR <- p.adjust(de$PValue, "BH")
#' plot_volcano(de)
#' plot_volcano(de, label = "selected", label.features = c("Gene1", "Gene2"))
#' plot_volcano(de, annotate.counts = "arrow", logfc.cutoff = 1)
#' \dontrun{
#' save_plot(plot_volcano(de), "volcano.pdf")
#' }
#'
#' @seealso \code{\link{standardize_de}}, \code{\link{save_plot}},
#'          \code{\link{export_source_data}}
#'
#' @import ggplot2
#' @export
#'
plot_volcano <- function(
    data,
    method = c("auto", "edgeR", "DESeq2", "limma", "seurat"),
    feature.col = NULL,
    logfc.col = NULL,
    p.col = NULL,
    padj.col = NULL,
    padj.cutoff = 0.05,
    logfc.cutoff = 0,
    label = c("top", "selected", "none"),
    label.n = 10,
    label.features = NULL,
    label.rank = c("padj", "logfc", "product"),
    label.size = NULL,
    italic.genes = TRUE,
    annotate.counts = c("corner", "arrow", "none"),
    symmetric = TRUE,
    cols = NULL,
    up.label = "Up-regulated",
    down.label = "Down-regulated",
    ns.label = "Not significant",
    pt.size = 1,
    alpha = 0.6,
    raster = NULL,
    raster.dpi = 1000,
    panel.width = 3,
    panel.height = panel.width,
    legend.position = NULL,
    bg = "transparent",
    seed = 42
) {
    method <- match.arg(method)
    label <- match.arg(label)
    label.rank <- match.arg(label.rank)
    annotate.counts <- match.arg(annotate.counts)
    if (!(is.numeric(padj.cutoff) && length(padj.cutoff) == 1 && padj.cutoff > 0 && padj.cutoff <= 1)) {
        stop("`padj.cutoff` must be a single number in (0, 1]", call. = FALSE)
    }
    if (!(is.numeric(logfc.cutoff) && length(logfc.cutoff) == 1 && logfc.cutoff >= 0)) {
        stop("`logfc.cutoff` must be a single number >= 0", call. = FALSE)
    }
    class_labels <- c(up = up.label, down = down.label, ns = ns.label)
    if (anyDuplicated(class_labels)) {
        stop("`up.label`, `down.label` and `ns.label` must be different", call. = FALSE)
    }
    cols <- volcano_colours(cols)
    if (is.null(legend.position)) {
        legend.position <- if (annotate.counts == "arrow") "none" else "right"
    }

    # Standard columns (feature, logfc, pvalue, padj, ...) from any DE tool
    df <- standardize_de(
        data,
        method = method,
        feature.col = feature.col,
        logfc.col = logfc.col,
        p.col = p.col,
        padj.col = padj.col
    )
    if (!is.null(df$comparison)) {
        stop("`plot_volcano()` draws one comparison; give a single DE table", call. = FALSE)
    }
    # A point needs both coordinates; NA padj (DESeq2 independent filtering)
    # or NA logFC (e.g. all-zero genes) cannot be placed on the plot
    drawable <- !is.na(df$logfc) & !is.na(df$padj)
    if (!all(drawable)) {
        message(sum(!drawable), " feature(s) without a log fold change or adjusted p-value are not drawn")
        df <- df[drawable, , drop = FALSE]
    }
    if (nrow(df) == 0) stop("No feature has both a log fold change and an adjusted p-value", call. = FALSE)

    # Direction of change: significance is decided by the adjusted p-value,
    # the same quantity as on the y-axis
    sig <- df$padj < padj.cutoff
    direction <- ifelse(sig & df$logfc > logfc.cutoff, "up",
                        ifelse(sig & df$logfc < -logfc.cutoff, "down", "ns"))
    df$direction <- factor(class_labels[direction], levels = class_labels)
    df$neg_log10_padj <- -log10(df$padj)
    df$labelled <- volcano_label_flags(df, direction, label, label.n, label.features, label.rank)
    n_up <- sum(direction == "up")
    n_down <- sum(direction == "down")

    # Not significant first, then down, then up: DE points are drawn last so
    # that the grey cloud never hides them
    df <- df[order(match(direction, c("ns", "down", "up"))), , drop = FALSE]
    rownames(df) <- NULL

    sz <- theme_sizes()
    lw <- pt_to_linewidth(sz$line)
    if (is.null(label.size)) label.size <- pt_to_size(sz$small)
    count_size <- pt_to_size(sz$base)

    # x-axis: symmetric around 0 so that up and down are comparable at a glance
    xlim <- if (symmetric) c(-1, 1) * max(abs(df$logfc)) else range(df$logfc)
    # y-axis: from 0, always reaching the significance line. The panel has a
    # fixed height, so space above the points (for the arrows) can be given in
    # inches and converted to data units
    headroom <- if (annotate.counts == "arrow") 3.2 * sz$small / 72 else 0
    ylim <- volcano_y_limits(
        y.max = max(df$neg_log10_padj, -log10(padj.cutoff)),
        panel.height = panel.height,
        headroom = headroom
    )

    layer_fun <- points_layer_fun(raster, nrow(df), 5000, raster.dpi)
    p <- ggplot(df, aes(x = .data$logfc, y = .data$neg_log10_padj)) +
        layer_fun(geom_point(aes(colour = .data$direction), size = pt.size, alpha = alpha,
                             shape = 16, stroke = 0)) +
        geom_hline(yintercept = -log10(padj.cutoff), linetype = "dashed", colour = "grey40",
                   linewidth = lw)
    if (logfc.cutoff > 0) {
        # Up to the top of the points only, so that the lines do not cross the
        # count arrows above them
        p <- p + annotate("segment", x = c(-1, 1) * logfc.cutoff, y = -Inf, yend = ylim$points.top,
                          linetype = "dashed", colour = "grey40", linewidth = lw)
    }
    p <- p +
        scale_colour_manual(values = stats::setNames(cols, class_labels), name = NULL) +
        scale_x_continuous(limits = xlim, expand = expansion(mult = 0.05)) +
        scale_y_continuous(limits = ylim$limits, expand = expansion(0)) +
        labs(x = quote(log[2] ~ "fold change"), y = quote(-log[10] ~ "adjusted" ~ italic(P))) +
        theme_bioToolkit("axis", bg = bg) +
        panel_size_theme(panel.width, panel.height) +
        theme(legend.position = legend.position) +
        # Legend keys at full opacity and a readable size, whatever the points
        guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))

    if (any(df$labelled)) {
        check_installed("ggrepel", reason = "to label features")
        p <- p + ggrepel::geom_text_repel(
            data = df[df$labelled, , drop = FALSE],
            aes(label = .data$feature),
            size = label.size,
            # Bold so that labels stay legible over the points and at print
            # sizes (user decision 2026-10-13); italic for gene names
            fontface = if (italic.genes) "bold.italic" else "bold",
            colour = "black",
            # Every requested label is drawn (ggrepel drops crowded labels by
            # default, which would silently lose selected features)
            max.overlaps = Inf,
            min.segment.length = 0,
            segment.colour = "grey30",
            segment.size = lw,
            box.padding = 0.25,
            # More time than ggrepel's default (0.5 s), which can stop before
            # crowded labels are separated
            max.time = 3,
            max.iter = 1e5,
            # Labels stay below the arrows, which sit above the points
            ylim = c(NA, ylim$points.top),
            seed = seed
        )
    }

    if (annotate.counts == "corner") {
        p <- p +
            annotate("text", x = Inf, y = -Inf, label = n_up, hjust = 1.2, vjust = -0.4,
                     size = count_size, colour = cols[["up"]]) +
            annotate("text", x = -Inf, y = -Inf, label = n_down, hjust = -0.2, vjust = -0.4,
                     size = count_size, colour = cols[["down"]])
    } else if (annotate.counts == "arrow") {
        p <- p + volcano_arrow_layers(
            xlim = xlim,
            y = ylim$arrow.y,
            labels = class_labels[c("down", "up")],
            counts = c(n_down, n_up),
            cols = cols[c("down", "up")],
            text.size = pt_to_size(sz$small),
            count.size = count_size,
            linewidth = pt_to_linewidth(sz$line * 1.5),
            head.length = sz$small * 0.5
        )
    }
    p
}


#' Colours of the three classes of a volcano plot
#'
#' @param cols NULL, colours named `up`, `down` and / or `ns` (the others keep
#'             their default), or three unnamed colours in that order
#'
#' @return A vector of three colours named `up`, `down`, `ns`
#'
#' @keywords internal
volcano_colours <- function(cols) {
    out <- c(up = "firebrick4", down = "midnightblue", ns = "grey75")
    if (is.null(cols)) return(out)
    if (!all(is_colour(cols))) stop("`cols` must be colours", call. = FALSE)
    if (is.null(names(cols))) {
        if (length(cols) != 3) {
            stop("Unnamed `cols` must be three colours: up, down, not significant", call. = FALSE)
        }
        names(cols) <- names(out)
    }
    unknown <- setdiff(names(cols), names(out))
    if (length(unknown)) {
        stop("`cols` must be named 'up', 'down' and / or 'ns', not: ", paste(unknown, collapse = ", "),
             call. = FALSE)
    }
    out[names(cols)] <- cols
    out
}


#' Which features of a volcano plot are labelled
#'
#' @param df Standardized DE table with `feature`, `logfc`, `padj`
#'
#' @param direction "up", "down" or "ns" per row of `df`
#'
#' @param label "top", "selected" or "none"
#'
#' @param label.n Number of top features per direction
#'
#' @param label.features Features to label, or NULL
#'
#' @param label.rank "padj", "logfc" or "product"
#'
#' @return A logical vector along the rows of `df`
#'
#' @keywords internal
volcano_label_flags <- function(df,
                                direction,
                                label,
                                label.n,
                                label.features,
                                label.rank) {
    flags <- rep(FALSE, nrow(df))
    if (label == "none") {
        if (!is.null(label.features)) message("`label.features` is ignored with `label = \"none\"`")
        return(flags)
    }
    if (label == "selected" && length(label.features) == 0) {
        stop("`label = \"selected\"` needs `label.features`", call. = FALSE)
    }

    if (label == "top") {
        if (!(is.numeric(label.n) && length(label.n) == 1 && label.n >= 0)) {
            stop("`label.n` must be a single number >= 0", call. = FALSE)
        }
        # Larger score = labelled first; order() breaks ties by the second key
        neg_log_p <- -log10(df$padj)
        abs_fc <- abs(df$logfc)
        for (d in c("up", "down")) {
            idx <- which(direction == d)
            ord <- switch(label.rank,
                padj = order(-neg_log_p[idx], -abs_fc[idx]),
                logfc = order(-abs_fc[idx], -neg_log_p[idx]),
                product = order(-(neg_log_p[idx] * abs_fc[idx]))
            )
            flags[idx[utils::head(ord, label.n)]] <- TRUE
        }
    }

    if (length(label.features)) {
        label.features <- as.character(label.features)
        missing <- setdiff(label.features, df$feature)
        if (length(missing)) {
            # Gene symbols differ in case between species (human CD3E, mouse
            # Cd3e); point to a case mismatch rather than failing silently
            near <- df$feature[tolower(df$feature) %in% tolower(missing)]
            warning("`label.features` not found among the drawn features: ",
                    paste(missing, collapse = ", "),
                    if (length(near)) paste0(" (matching only in case: ", paste(unique(near), collapse = ", "), ")"),
                    call. = FALSE)
        }
        flags <- flags | df$feature %in% label.features
    }
    flags
}


#' y-axis limits of a volcano plot
#'
#' The panel has a fixed height, so a given space in inches above the points
#' can be converted to data units.
#'
#' @param y.max Highest -log10 adjusted p-value to show
#'
#' @param panel.height Panel height in inches
#'
#' @param headroom Space in inches to keep free above the points
#'
#' @return A list with the axis `limits`, the top of the points region
#'         (`points.top`, the upper limit for labels) and the height of the
#'         arrows (`arrow.y`, the middle of the headroom)
#'
#' @keywords internal
volcano_y_limits <- function(y.max,
                             panel.height,
                             headroom) {
    # Avoid a zero-height axis when nothing approaches significance
    y.max <- max(y.max, 1e-3)
    # Margins of 3% below (points at y = 0 are not cut by the axis) and 5%
    # above the points, as fractions of the panel
    bottom <- 0.03 * panel.height
    top <- 0.05 * panel.height
    if (headroom >= panel.height - bottom - top) {
        stop("`panel.height` is too small for the count arrows", call. = FALSE)
    }
    per_inch <- y.max / (panel.height - bottom - top - headroom)
    points_top <- y.max + top * per_inch
    list(
        limits = c(-bottom * per_inch, points_top + headroom * per_inch),
        points.top = points_top,
        arrow.y = points_top + headroom / 2 * per_inch
    )
}


#' Arrows naming the directions of change, with counts, for a volcano plot
#'
#' Two arrows at height `y`, pointing left (down-regulated) and right
#' (up-regulated) from the centre of the x range, each with its label above
#' and its count below.
#'
#' @param xlim x range the arrows are laid out in (the axis limits)
#'
#' @param y Height of the arrows (data units)
#'
#' @param labels Labels of the down and up arrows
#'
#' @param counts Counts of down- and up-regulated features
#'
#' @param cols Colours of the down and up arrows
#'
#' @param text.size Size of the labels (mm)
#'
#' @param count.size Size of the counts (mm)
#'
#' @param linewidth Line width of the arrows (ggplot2 units)
#'
#' @param head.length Length of the arrow heads in points
#'
#' @return A list of ggplot2 layers
#'
#' @keywords internal
volcano_arrow_layers <- function(xlim,
                                 y,
                                 labels,
                                 counts,
                                 cols,
                                 text.size,
                                 count.size,
                                 linewidth,
                                 head.length) {
    w <- diff(xlim)
    centre <- mean(xlim)
    # Each arrow covers 40% of the width, from 5% off-centre outwards, so both
    # stay inside the axis limits (and are not dropped by the scale)
    inner <- 0.05 * w
    outer <- 0.45 * w
    x_start <- centre + c(-inner, inner)
    x_end <- centre + c(-outer, outer)
    x_mid <- (x_start + x_end) / 2
    head <- grid::arrow(length = grid::unit(head.length, "pt"), type = "closed")
    list(
        annotate("segment", x = x_start, xend = x_end, y = y, yend = y, colour = cols,
                 linewidth = linewidth, arrow = head, arrow.fill = cols),
        annotate("text", x = x_mid, y = y, label = labels, vjust = -0.5, size = text.size, colour = cols),
        annotate("text", x = x_mid, y = y, label = counts, vjust = 1.5, size = count.size, colour = cols)
    )
}
