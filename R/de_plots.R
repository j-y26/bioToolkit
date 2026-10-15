# Summaries of differential expression across comparisons (roadmap items
# G2 and G4)
#
# Design
# ------
# Both plots read many comparisons at once through `standardize_de()` (F5):
# a named list of DE tables, or one long table with a comparison column
# (e.g. FindAllMarkers() output with `comparison.col = "cluster"`), with an
# optional second level of `groups` (e.g. cell type above contrasts).
#
# * `plot_de_dots()` (G2) shows *selected* genes: a grid of genes x
#   comparisons on the dot-grid engine (G0), colour = log2 fold change
#   (diverging, centred on 0, each side to its own extreme), size = -log10
#   adjusted p-value. Comparisons are marked by a coloured dot, as the groups
#   of DotPlot_bioToolkit(); `groups` and gene groups by bars (R/group_bars.R). Significance
#   needs both cutoffs (padj and |log2FC|), so it is not readable from the
#   dot size alone: significant results are drawn opaque with a dark ring,
#   the others faded (default, "fade"), keeping their fold-change colour so
#   that no information is lost. One titled legend ("Significance") spells
#   out both cutoffs. "grey" and "hide" remain available. Replaces
#   `de_dot_plot()`.
# * `plot_multi_volcano()` (G4) shows *all* DE genes: one column of genes per
#   comparison, log2 fold change on the y-axis. No capping of extreme values
#   (decided 2026-10-13). Revision of 2026-10-14 (user): the coloured band at
#   0 spanning +-`logfc.cutoff` is gone (huge with a high cutoff, invisible at
#   0); each comparison is marked by a strip of fixed thickness above its
#   column (R/group_bars.R) with its name beyond, the cutoffs are dashed
#   lines, non-significant genes are drawn in grey by default (`show.ns`),
#   `groups` get bars like split.by, and `flip = TRUE` lists comparisons as
#   rows (names on the left, down-regulated genes left, up right).


#' @title Dot plot of differential expression across comparisons
#'
#' @description Shows selected genes across several DE comparisons: each dot's
#'              colour is the log2 fold change and its size the -log10
#'              adjusted p-value; non-significant results are faded. Reads
#'              edgeR, DESeq2, limma and Seurat results (see
#'              \code{\link{standardize_de}}). Every row and column has a fixed
#'              size, so the plot is sized the same way in every figure; save
#'              it with \code{\link{save_plot}}.
#'
#' @param de DE results of several comparisons: a named list of DE tables
#'           (names = comparisons), or one table with a comparison column
#'           (`comparison.col`), e.g. \code{Seurat::FindAllMarkers()} output
#'
#' @param features Genes to show: a character vector, or a named list (gene
#'                 groups, drawn as labelled groups of columns). If NULL
#'                 (default), the `n` most significant genes of each
#'                 comparison (by adjusted p-value, then absolute log2 fold
#'                 change), combined. Matched exactly (case matters)
#'
#' @param n Number of top genes per comparison when `features` is NULL,
#'          default is 5
#'
#' @param method,feature.col,logfc.col,p.col,padj.col,comparison.col,groups
#'        How the DE tables are read, see \code{\link{standardize_de}}.
#'        `groups` (a group per comparison) draws comparisons as labelled
#'        groups of rows
#'
#' @param padj.cutoff,logfc.cutoff A result is significant when its adjusted
#'                                 p-value is below `padj.cutoff` (default
#'                                 0.05) and its absolute log2 fold change
#'                                 above `logfc.cutoff` (default 0)
#'
#' @param significance How significance is shown: "fade" (default;
#'                     significant dots opaque with a dark ring, the others
#'                     faded, all keeping their fold-change colour), "grey"
#'                     (non-significant dots grey) or "hide" (non-significant
#'                     dots not drawn)
#'
#' @param size.max -log10 adjusted p-value drawn as the largest dot, default
#'                 is 20 (padj = 1e-20); more significant results are drawn at
#'                 that size and the legend says so
#'
#' @param logfc.limits Limits of the colour scale. If NULL (default), from
#'                     the most negative to the most positive log2 fold
#'                     change of the significant dots, centred on 0: each
#'                     side reaches its full colour at its own extreme. Give
#'                     two numbers to share one scale between figures; values
#'                     beyond get the colour of the end
#'
#' @param cluster.features,cluster.comparisons Whether to order genes (within
#'        their groups) / comparisons (within their groups) by hierarchical
#'        clustering of their log2 fold changes, default is FALSE
#'
#' @param cols Colours of the diverging log2 fold change scale: NULL
#'             (default, "blue-magenta"), a palette name, or colours from low
#'             to high, e.g. `c("#2166AC", "white", "#B2182B")`
#'
#' @param comparison.marker How each comparison is marked with its colour
#'                          next to its name, as the groups of
#'                          \code{\link{DotPlot_bioToolkit}}: "dot" (default;
#'                          a coloured dot), "strip" (a colour strip) or
#'                          "none"
#'
#' @param comparison.cols Colours of the comparison markers: NULL (default,
#'                        "scanpy-20"), a qualitative palette name, or colours
#'                        named by comparison (e.g. the cell type colours of
#'                        an embedding)
#'
#' @param groups.cols Colours of the bars marking the `groups` of comparisons:
#'                    NULL (default) for black bars, a qualitative palette
#'                    name, or colours named by group
#'
#' @param groups.side Side of the `groups`: "opposite" (default; opposite the
#'                    comparison names: on the right, or on top with `flip =
#'                    TRUE`) or "labels" (beyond the comparison names)
#'
#' @param feature.group.cols Colours of the bars marking the groups of genes
#'                           (a named list of `features`): NULL (default) for
#'                           black bars, a qualitative palette name, or colours
#'                           named by gene group
#'
#' @param feature.group.side Side of the gene groups: "opposite" (default; on
#'                           top, or on the right with `flip = TRUE`) or
#'                           "labels" (beyond the gene names)
#'
#' @param group.bar.width Thickness of the bars marking gene groups and
#'                        `groups`, in points, default is 3.5
#'
#' @param flip Whether to draw genes as rows and comparisons as columns,
#'             default is FALSE
#'
#' @param italic.genes Whether gene names are set in italics, default is TRUE
#'
#' @param cell.width,cell.height Width of a column and height of a row in
#'                               inches, default is 0.18
#'
#' @param legend.position Position of the legends, default is "right"
#'
#' @param bg Background colour, default is "transparent"
#'
#' @details Comparisons play the part of `group.by` in the single-cell plots
#'          (the units compared, marked by colour) and `groups` that of
#'          `split.by` (higher-level blocks, marked by a bar and a name).
#'          Genes not tested in a comparison have no dot. The source data
#'          (\code{\link{export_source_data}}) hold one row per gene and
#'          comparison with the standardized DE columns plus `significant`.
#'
#' @return A ggplot object
#'
#' @examples
#' set.seed(1)
#' genes <- paste0("Gene", 1:200)
#' make_de <- function() {
#'     d <- data.frame(logFC = rnorm(200, sd = 1.5), logCPM = 5, PValue = runif(200)^5,
#'                     row.names = genes)
#'     d$FDR <- p.adjust(d$PValue, "BH")
#'     d
#' }
#' de <- list("KO vs WT" = make_de(), "Het vs WT" = make_de(), "KO vs Het" = make_de())
#' plot_de_dots(de, n = 4)
#' plot_de_dots(de, features = list(Up = c("Gene1", "Gene2"), Down = c("Gene3", "Gene4")))
#'
#' @seealso \code{\link{standardize_de}}, \code{\link{plot_multi_volcano}},
#'          \code{\link{plot_volcano}}
#'
#' @import ggplot2
#' @export
#'
plot_de_dots <- function(
    de,
    features = NULL,
    n = 5,
    method = c("auto", "edgeR", "DESeq2", "limma", "seurat"),
    feature.col = NULL,
    logfc.col = NULL,
    p.col = NULL,
    padj.col = NULL,
    comparison.col = NULL,
    groups = NULL,
    padj.cutoff = 0.05,
    logfc.cutoff = 0,
    significance = c("fade", "grey", "hide"),
    size.max = 20,
    logfc.limits = NULL,
    cluster.features = FALSE,
    cluster.comparisons = FALSE,
    cols = NULL,
    comparison.marker = c("dot", "strip", "none"),
    comparison.cols = NULL,
    groups.cols = NULL,
    groups.side = c("opposite", "labels"),
    feature.group.cols = NULL,
    feature.group.side = c("opposite", "labels"),
    group.bar.width = 3.5,
    flip = FALSE,
    italic.genes = TRUE,
    cell.width = 0.18,
    cell.height = 0.18,
    legend.position = c("right", "left", "top", "bottom", "none"),
    bg = "transparent"
) {
    method <- match.arg(method)
    significance <- match.arg(significance)
    comparison.marker <- match.arg(comparison.marker)
    groups.side <- match.arg(groups.side)
    feature.group.side <- match.arg(feature.group.side)
    legend.position <- match.arg(legend.position)
    check_de_cutoffs(padj.cutoff, logfc.cutoff)
    if (!(is.numeric(size.max) && length(size.max) == 1 && size.max > 0)) {
        stop("`size.max` must be a single positive number", call. = FALSE)
    }

    res <- standardize_de(de, method = method, feature.col = feature.col, logfc.col = logfc.col,
                          p.col = p.col, padj.col = padj.col, comparison.col = comparison.col,
                          groups = groups)
    if (is.null(res$comparison)) {
        res$comparison <- factor("Comparison")
        message("A single comparison: give a named list of tables, or `comparison.col`, for several")
    }
    res$significant <- !is.na(res$padj) & res$padj < padj.cutoff &
        !is.na(res$logfc) & abs(res$logfc) > logfc.cutoff

    # Genes: given, or the top `n` of each comparison
    if (is.null(features)) {
        sig <- res[res$significant, , drop = FALSE]
        if (nrow(sig) == 0) stop("No significant results to pick top genes from; give `features`", call. = FALSE)
        sig <- sig[order(sig$padj, -abs(sig$logfc)), , drop = FALSE]
        features <- unique(unlist(lapply(split(sig$feature, sig$comparison), utils::head, n)))
    }
    ft <- dot_feature_table(features)
    missing <- setdiff(unique(ft$feature), res$feature)
    if (length(missing)) {
        near <- unique(res$feature[tolower(res$feature) %in% tolower(missing)])
        warning(length(missing), " feature(s) not in the DE results are left out: ",
                paste(missing, collapse = ", "),
                if (length(near)) paste0(" (matching only in case: ", paste(near, collapse = ", "), ")"),
                call. = FALSE)
        ft <- ft[!(ft$feature %in% missing), , drop = FALSE]
    }
    if (nrow(ft) == 0) stop("None of `features` are in the DE results", call. = FALSE)

    # One row per gene (as listed, a gene in two groups twice) and comparison
    cmp <- levels(res$comparison)
    idx <- expand.grid(j = seq_len(nrow(ft)), i = seq_along(cmp))
    df <- data.frame(id = ft$id[idx$j], feature = ft$feature[idx$j], comparison = cmp[idx$i],
                     stringsAsFactors = FALSE)
    if (!is.null(ft$feature_group)) df$feature_group <- ft$feature_group[idx$j]
    key <- paste(res$comparison, res$feature, sep = "\r")
    hit <- match(paste(df$comparison, df$feature, sep = "\r"), key)
    df$logfc <- res$logfc[hit]
    df$padj <- res$padj[hit]
    df$significant <- res$significant[hit] %in% TRUE
    if (!is.null(res$group)) df$group <- res$group[match(df$comparison, as.character(res$comparison))]
    df$neg_log10_padj <- -log10(df$padj)
    df$comparison <- factor(df$comparison, levels = cmp)

    # Orders
    feat_lv <- unique(ft$id)
    cmp_lv <- cmp
    if (cluster.features || cluster.comparisons) {
        m <- tapply(df$logfc, list(df$id, as.character(df$comparison)), function(v) v[1])
        if (cluster.features) feat_lv <- cluster_order(m[feat_lv, cmp, drop = FALSE], ft$feature_group)
        if (cluster.comparisons) {
            # Comparisons stay within their groups
            grp <- if (!is.null(df$group)) df$group[match(cmp, as.character(df$comparison))]
            cmp_lv <- cluster_order(t(m[feat_lv, cmp, drop = FALSE]), grp)
        }
    }

    # Dots: size from the adjusted p-value; NA when not tested (or hidden)
    drawn <- !is.na(df$logfc) & !is.na(df$padj)
    if (significance == "hide") drawn <- drawn & df$significant
    size <- ifelse(drawn, df$neg_log10_padj, NA)
    fill <- if (significance == "grey") ifelse(df$significant, df$logfc, NA) else df$logfc
    if (is.null(logfc.limits)) {
        # From the most negative to the most positive fold change of the
        # significant dots (all dots when none is), centred on 0: each side
        # gets its full colour at its own extreme (diverging_limits())
        v <- df$logfc[drawn & df$significant]
        if (length(v) == 0) v <- df$logfc[drawn]
        logfc.limits <- diverging_limits(v, center = 0)
    }
    size_breaks <- scales::breaks_extended(4)(c(0, size.max))
    size_breaks <- size_breaks[size_breaks > 0 & size_breaks <= size.max]
    size_labels <- format(size_breaks, trim = TRUE)
    if (any(size > size.max, na.rm = TRUE)) {
        if (utils::tail(size_breaks, 1) != size.max) {
            size_breaks <- c(size_breaks, size.max)
            size_labels <- c(size_labels, format(size.max, trim = TRUE))
        }
        size_labels[length(size_labels)] <- paste0("\u2265", size_labels[length(size_labels)])
    }

    feat <- factor(df$id, levels = feat_lv)
    comp <- factor(as.character(df$comparison), levels = cmp_lv)
    feat_group <- if (!is.null(ft$feature_group)) factor(df$feature_group, levels = levels(ft$feature_group))
    comp_group <- if (!is.null(df$group)) factor(as.character(df$group), levels = levels(res$group))
    eng <- data.frame(size = size, fill = fill)
    if (!flip) {
        eng$x <- feat
        eng$y <- comp
        eng$x_group <- feat_group
        eng$y_group <- comp_group
    } else {
        eng$x <- comp
        eng$y <- feat
        eng$x_group <- comp_group
        eng$y_group <- feat_group
    }
    if (significance == "fade") eng$outline <- df$significant
    # Comparisons are marked by their colour (as the groups of DotPlot)
    has_marker <- comparison.marker != "none"
    if (has_marker) {
        cmp_cols <- group_colours(comparison.cols, factor(cmp, levels = cmp), arg = "comparison.cols")$cols
        eng$strip <- unname(cmp_cols[as.character(df$comparison)])
    }
    feat_labels <- stats::setNames(ft$feature, ft$id)
    cmp_labels <- stats::setNames(cmp, cmp)
    italic <- stats::setNames(rep(italic.genes, length(feat_labels)), names(feat_labels))
    # Bar colours of the gene groups and comparison groups (black by default)
    feat_bar_cols <- if (!is.null(feat_group)) {
        group_bar_colours(feature.group.cols, levels(droplevels(feat_group)), "feature.group.cols")
    }
    comp_bar_cols <- if (!is.null(comp_group)) {
        group_bar_colours(groups.cols, levels(droplevels(comp_group)), "groups.cols")
    }

    p <- dot_grid_plot(
        eng,
        x.labels = if (!flip) feat_labels else cmp_labels,
        y.labels = if (!flip) cmp_labels else feat_labels,
        x.italic = if (!flip) italic else FALSE,
        y.italic = if (!flip) FALSE else italic,
        size.limits = c(0, size.max),
        size.breaks = size_breaks,
        size.labels = size_labels,
        size.title = quote(-log[10] ~ "padj"),
        fill.colours = resolve_continuous_palette(cols, diverging = TRUE),
        fill.limits = logfc.limits,
        fill.center = 0,
        fill.title = quote(log[2] ~ "fold change"),
        na.value = "grey80",
        outline.labels = significance_labels(padj.cutoff, logfc.cutoff),
        outline.title = "Significance",
        strip.axis = if (has_marker) (if (!flip) "y" else "x"),
        strip.style = if (comparison.marker == "dot") "dot" else "strip",
        cell.width = cell.width,
        cell.height = cell.height,
        x.group.side = if (!flip) feature.group.side else groups.side,
        y.group.side = if (!flip) groups.side else feature.group.side,
        x.group.cols = if (!flip) feat_bar_cols else comp_bar_cols,
        y.group.cols = if (!flip) comp_bar_cols else feat_bar_cols,
        bar.width = group.bar.width,
        legend.position = legend.position,
        bg = bg
    )
    src <- df[, intersect(c("group", "comparison", "feature_group", "feature", "logfc", "padj",
                            "neg_log10_padj", "significant"), names(df))]
    attr(p, "source_data") <- src
    p
}


#' Legend labels of significant and other results
#'
#' Both cutoffs are spelled out, since a result can have a small adjusted
#' p-value and still not be significant for its fold change.
#'
#' @param padj.cutoff Adjusted p-value cutoff
#'
#' @param logfc.cutoff Absolute log2 fold change cutoff (0: none)
#'
#' @return Two labels: significant, not significant
#'
#' @keywords internal
significance_labels <- function(padj.cutoff,
                                logfc.cutoff = 0) {
    c(paste0("padj < ", format(padj.cutoff),
             if (logfc.cutoff > 0) paste0("\nand |log2FC| > ", format(logfc.cutoff))),
      "Not significant")
}


#' Check the significance cutoffs of DE plots
#'
#' @param padj.cutoff Adjusted p-value cutoff
#'
#' @param logfc.cutoff Absolute log2 fold change cutoff
#'
#' @keywords internal
check_de_cutoffs <- function(padj.cutoff,
                             logfc.cutoff) {
    if (!(is.numeric(padj.cutoff) && length(padj.cutoff) == 1 && padj.cutoff > 0 && padj.cutoff <= 1)) {
        stop("`padj.cutoff` must be a single number in (0, 1]", call. = FALSE)
    }
    if (!(is.numeric(logfc.cutoff) && length(logfc.cutoff) == 1 && logfc.cutoff >= 0)) {
        stop("`logfc.cutoff` must be a single number >= 0", call. = FALSE)
    }
    invisible(TRUE)
}


#' @title Volcano-style summary of many DE comparisons
#'
#' @description Summarizes genome-wide DE results of many comparisons in one
#'              figure (e.g. condition B vs A within each cell type): one
#'              column of genes per comparison, at their log2 fold change
#'              (non-significant genes in grey, the cutoffs as dashed lines),
#'              each comparison marked by a coloured strip and its name on
#'              top, the numbers of up- and down-regulated genes above and
#'              below each column, and the top genes labelled. `flip = TRUE`
#'              lists the comparisons as rows. An improved version of the plot made
#'              popular by scRNAtoolVis' \code{jjVolcano()}: colour carries
#'              the significance (or how many comparisons share a gene)
#'              rather than repeating the direction, points are spread by
#'              their density so that each column shows the shape of its
#'              distribution, and no value is capped.
#'
#' @param de DE results of several comparisons: a named list of DE tables
#'           (names = comparisons), or one table with a comparison column
#'           (`comparison.col`), e.g. \code{Seurat::FindAllMarkers()} output
#'           with `comparison.col = "cluster"`. For single-cell data,
#'           pseudobulk DE results are recommended: per-cell tests (e.g.
#'           Wilcoxon) treat cells as replicates and inflate the numbers of DE
#'           genes, more so in large groups
#'
#' @param method,feature.col,logfc.col,p.col,padj.col,comparison.col,groups
#'        How the DE tables are read, see \code{\link{standardize_de}}.
#'        `groups` (a group per comparison) draws the columns of a group
#'        together under its label
#'
#' @param padj.cutoff,logfc.cutoff A gene is DE when its adjusted p-value is
#'                                 below `padj.cutoff` (default 0.05) and its
#'                                 absolute log2 fold change above
#'                                 `logfc.cutoff` (default 1)
#'
#' @param show.ns Whether to draw non-significant genes too, in grey under
#'                the DE genes as in \code{\link{plot_volcano}}, default is
#'                TRUE
#'
#' @param layout How points spread across their column: "sina" (default;
#'               as wide as the density of genes at that fold change, so the
#'               shape of the distribution shows), "pvalue" (x within the
#'               column = -log10 adjusted p-value, from the left, making each
#'               column a sideways volcano) or "jitter" (uniform, the classic
#'               look)
#'
#' @param color.by What the colour of a point shows: "padj" (default; -log10
#'                 adjusted p-value on two ramps, warm for up and cool for
#'                 down, sharing one set of limits; the colour saturates at
#'                 the 95th percentile of the DE genes' values, so that a few
#'                 p-values near 0 do not leave the others pale, while the
#'                 points themselves are never capped), "direction" (one colour
#'                 for up and one for down) or "sharing" (in how many of the
#'                 comparisons the gene is DE: comparison-specific versus
#'                 shared responses)
#'
#' @param cols Colours of the points: for "padj", a diverging palette name or
#'             colours from low (down) to high (up); for "direction", colours
#'             named `up` / `down` (or two colours, up then down); for
#'             "sharing", a sequential palette name or colours from low to
#'             high. NULL (default) uses "blue-magenta" without its pale
#'             centre, `c(up = "firebrick4", down = "midnightblue")` and
#'             "plasma_r"
#'
#' @param comparison.cols Colours of the comparison strips (a bar along
#'                        each column, or row when flipped): NULL (default;
#'                        the default group palette, "scanpy-20"), a
#'                        qualitative palette name, or colours named by
#'                        comparison (e.g. the cell type colours of an
#'                        embedding)
#'
#' @param comparison.labels Display names of the comparisons, named by
#'                          comparison, written beyond the strips (turned 45
#'                          degrees when wider than a column)
#'
#' @param groups.cols Colours of the bars marking the `groups` of comparisons:
#'                    NULL (default) for black bars, a qualitative palette
#'                    name, or colours named by group
#'
#' @param groups.side Side of the `groups`: "opposite" (default; opposite the
#'                    comparison names: at the bottom, or on the right when
#'                    flipped) or "labels" (beyond the comparison names)
#'
#' @param group.bar.width Thickness of the comparison strips and group bars in
#'                        points, default is 3.5
#'
#' @param flip Whether to list the comparisons as rows, with log2 fold change
#'             on the x-axis (down-regulated genes on the left, up on the
#'             right) and the comparison names on the left, default is FALSE
#'             (comparisons as columns, names on top)
#'
#' @param label Which genes to label: "top" (default; the `label.n` most
#'              significant genes per column and direction, plus
#'              `label.features`), "selected" (only `label.features`) or
#'              "none"
#'
#' @param label.n Number of genes labelled per column and direction, default
#'                is 3
#'
#' @param label.features Genes always labelled (matched exactly)
#'
#' @param label.rank How the top genes are ranked: "padj" (default), "logfc"
#'                   or "product", as in \code{\link{plot_volcano}}
#'
#' @param show.counts Whether to write the numbers of up- and down-regulated
#'                    genes above and below each column, default is TRUE
#'
#' @param order Order of the comparisons: "input" (default; as given),
#'              "n_de" (most DE genes first) or "cluster" (similar log2 fold
#'              change profiles together). Within groups when `groups` is
#'              given
#'
#' @param italic.genes Whether labels are gene names, in italics, default is
#'                     TRUE
#'
#' @param pt.size Point size, default is 0.6
#'
#' @param alpha Point opacity, default is 0.8
#'
#' @param raster Whether to rasterize the points (labels, strips and text stay
#'               vector): NULL (default) rasterizes above 5,000 points, or
#'               TRUE / FALSE
#'
#' @param raster.dpi Resolution of rasterized points, default is 1000
#'
#' @param cell.width Width of one column in inches, default is 0.6
#'
#' @param cell.height Height of one row in inches when flipped, default is
#'                    `cell.width`
#'
#' @param panel.height Height of the panel (the log2 fold change axis) in
#'                     inches, default is 3.5
#'
#' @param panel.width Width of the panel (the log2 fold change axis) in inches
#'                    when flipped, default is `panel.height`
#'
#' @param legend.position Position of the legend, default is "right"
#'
#' @param bg Background colour, default is "transparent"
#'
#' @param seed Random seed for the spread of points and label placement,
#'             default is 42
#'
#' @param ... `band.cols` and `band.labels`, the names of `comparison.cols`
#'            and `comparison.labels` before 0.11.0, accepted with a
#'            deprecation warning
#'
#' @details Comparisons play the part of `group.by` in the single-cell plots
#'          (the units compared, marked by colour) and `groups` that of
#'          `split.by` (higher-level blocks, marked by a bar and a name).
#'          The source data (\code{\link{export_source_data}}) hold one row
#'          per drawn gene with its comparison, x position, log2 fold change,
#'          adjusted p-value, direction, the number of comparisons it is DE in
#'          (`n_shared`) and whether it is labelled.
#'
#' @return A ggplot object
#'
#' @examples
#' set.seed(1)
#' make_de <- function(shift) {
#'     d <- data.frame(logFC = rnorm(3000, sd = 1.2) + shift, logCPM = 5,
#'                     PValue = runif(3000)^6, row.names = paste0("Gene", 1:3000))
#'     d$FDR <- p.adjust(d$PValue, "BH")
#'     d
#' }
#' de <- list("T cells" = make_de(0), "B cells" = make_de(0.3), "Monocytes" = make_de(-0.2))
#' plot_multi_volcano(de)
#' plot_multi_volcano(de, layout = "pvalue", color.by = "direction")
#'
#' @seealso \code{\link{plot_de_dots}}, \code{\link{plot_volcano}},
#'          \code{\link{standardize_de}}
#'
#' @import ggplot2
#' @export
#'
plot_multi_volcano <- function(
    de,
    method = c("auto", "edgeR", "DESeq2", "limma", "seurat"),
    feature.col = NULL,
    logfc.col = NULL,
    p.col = NULL,
    padj.col = NULL,
    comparison.col = NULL,
    groups = NULL,
    padj.cutoff = 0.05,
    logfc.cutoff = 1,
    show.ns = TRUE,
    layout = c("sina", "pvalue", "jitter"),
    color.by = c("padj", "direction", "sharing"),
    cols = NULL,
    comparison.cols = NULL,
    comparison.labels = NULL,
    groups.cols = NULL,
    groups.side = c("opposite", "labels"),
    group.bar.width = 3.5,
    flip = FALSE,
    label = c("top", "selected", "none"),
    label.n = 3,
    label.features = NULL,
    label.rank = c("padj", "logfc", "product"),
    show.counts = TRUE,
    order = c("input", "n_de", "cluster"),
    italic.genes = TRUE,
    pt.size = 0.6,
    alpha = 0.8,
    raster = NULL,
    raster.dpi = 1000,
    cell.width = 0.6,
    cell.height = cell.width,
    panel.height = 3.5,
    panel.width = panel.height,
    legend.position = "right",
    bg = "transparent",
    seed = 42,
    ...
) {
    # `band.cols` / `band.labels` (before 0.11.0, when comparisons were
    # coloured bands) are the new `comparison.cols` / `comparison.labels`
    list2env(renamed_args(list(...), c(band.cols = "comparison.cols", band.labels = "comparison.labels"),
                          "plot_multi_volcano", names(match.call())), environment())
    method <- match.arg(method)
    layout <- match.arg(layout)
    color.by <- match.arg(color.by)
    label <- match.arg(label)
    label.rank <- match.arg(label.rank)
    order <- match.arg(order)
    groups.side <- match.arg(groups.side)
    check_bar_width(group.bar.width)
    check_de_cutoffs(padj.cutoff, logfc.cutoff)
    sz <- theme_sizes()
    lw <- pt_to_linewidth(sz$line)

    res <- standardize_de(de, method = method, feature.col = feature.col, logfc.col = logfc.col,
                          p.col = p.col, padj.col = padj.col, comparison.col = comparison.col,
                          groups = groups)
    if (is.null(res$comparison)) {
        stop("`plot_multi_volcano()` summarizes several comparisons; give a named list of tables ",
             "or `comparison.col` (use plot_volcano() for one)", call. = FALSE)
    }
    res <- res[!is.na(res$logfc) & !is.na(res$padj), , drop = FALSE]
    sig <- res$padj < padj.cutoff
    res$direction <- factor(ifelse(sig & res$logfc > logfc.cutoff, "up",
                                   ifelse(sig & res$logfc < -logfc.cutoff, "down", "ns")),
                            levels = c("up", "down", "ns"))
    res$neg_log10_padj <- -log10(res$padj)
    is_de <- res$direction != "ns"
    # In how many comparisons each gene is DE
    shared <- table(unique(res[is_de, c("feature", "comparison")])$feature)
    res$n_shared <- as.integer(shared[res$feature])
    res$n_shared[is.na(res$n_shared)] <- 0L

    # Order of comparisons, within groups
    cmp <- levels(res$comparison)
    grp_of <- if (!is.null(res$group)) res$group[match(cmp, as.character(res$comparison))]
    if (order == "n_de") {
        n_de <- tapply(is_de, res$comparison, sum)[cmp]
        cmp <- cmp[if (is.null(grp_of)) base::order(-n_de) else base::order(grp_of, -n_de)]
    } else if (order == "cluster") {
        genes <- unique(res$feature[is_de])
        m <- tapply(res$logfc, list(as.character(res$comparison), res$feature), function(v) v[1])
        m <- m[cmp, intersect(genes, colnames(m)), drop = FALSE]
        cmp <- cluster_order(m, grp_of)
    }
    grp_of <- if (!is.null(res$group)) res$group[match(cmp, as.character(res$comparison))]
    # The comparisons of a group sit together (in the order above), so that
    # each group's panel spans consecutive positions
    if (!is.null(grp_of)) {
        cmp <- cmp[base::order(grp_of)]
        grp_of <- grp_of[base::order(grp_of)]
    }

    df <- res[is_de | show.ns, , drop = FALSE]
    if (!any(df$direction != "ns")) warning("No DE genes at these cutoffs", call. = FALSE)
    df$comparison <- factor(as.character(df$comparison), levels = cmp)
    # Position of each comparison, numbered across all groups so that one
    # axis names them in every panel of a grouped plot (each group's panel
    # spans its own positions). Rows of a flipped plot list the first
    # comparison at the top
    pos <- stats::setNames(if (flip) rev(seq_along(cmp)) else seq_along(cmp), cmp)
    half <- 0.4
    df$pos <- pos[as.character(df$comparison)] +
        with_seed(seed, multi_volcano_offsets(df, layout, half))

    # Value (log2 fold change) limits: all values (no capping), with room in
    # inches at both ends for the counts
    v_rng <- range(c(df$logfc, -logfc.cutoff, logfc.cutoff, 0))
    room <- if (show.counts) (if (flip) 0.3 else 0.22) else 0.05
    value_len <- if (flip) panel.width else panel.height
    unit <- diff(v_rng) / (value_len - 2 * room)
    vlim <- v_rng + c(-1, 1) * room * unit

    # Comparisons: display names and strip colours
    cmp_text <- stats::setNames(cmp, cmp)
    if (!is.null(comparison.labels)) {
        known <- intersect(names(comparison.labels), cmp)
        cmp_text[known] <- comparison.labels[known]
    }
    cmp_cols <- group_colours(comparison.cols, factor(cmp, levels = cmp), arg = "comparison.cols")$cols
    # Names above the columns are turned 45 degrees when one is wider than
    # its column (about 0.55 em per character)
    upright <- flip || all(nchar(cmp_text) * 0.55 * sz$small / 72 <= 0.92 * cell.width)

    # Labels: top genes per comparison and direction, plus selected genes
    df$labelled <- FALSE
    if (label == "top") {
        for (cm in cmp) {
            idx <- which(df$comparison == cm)
            dir <- as.character(df$direction[idx])
            df$labelled[idx] <- volcano_label_flags(df[idx, , drop = FALSE], dir, "top", label.n, NULL, label.rank)
        }
    }
    if (label != "none" && length(label.features)) {
        df$labelled <- df$labelled | volcano_label_flags(df, as.character(df$direction), "selected", 0,
                                                         label.features, label.rank)
    } else if (label == "selected") {
        stop("`label = \"selected\"` needs `label.features`", call. = FALSE)
    }

    # Point colours
    up_down <- volcano_colours(if (color.by == "direction") cols else NULL)
    colourbar <- guide_colourbar(theme = theme(
        legend.key.height = unit(1.2, "in"),
        legend.key.width = unit(0.12, "in"),
        legend.frame = element_rect(colour = "black", linewidth = 0.3),
        legend.ticks = element_line(colour = "black", linewidth = 0.3)
    ))
    if (color.by == "padj") {
        df$colour_value <- ifelse(df$direction == "up", 1, -1) * df$neg_log10_padj
        # The colour saturates at the 95th percentile of the DE genes' -log10
        # padj: a few p-values near 0 (common with per-cell tests) would
        # otherwise leave every other gene pale. Positions are not capped.
        # Each end is trimmed to the most extreme gene of its direction, so
        # that the side with fewer or weaker genes is not faded
        # (diverging_limits(); 2026-10-14)
        v <- df$colour_value[df$direction != "ns"]
        m <- max(if (length(v)) stats::quantile(abs(v), 0.95, names = FALSE) else 0, -log10(padj.cutoff) * 2)
        lim <- diverging_limits(v, -m, m)
        sat <- c(any(v < lim[1]), any(v > lim[2]))
        pal <- if (is.null(cols)) bt_palette("blue-magenta")[c(1:4, 8:11)] else resolve_continuous_palette(cols, TRUE)
        # Breaks on each side up to its own end, which is labelled; a side
        # gets breaks in proportion to its share of the colour bar, so that
        # the labels of a short side do not crowd
        side_breaks <- function(end) {
            if (end <= 0) return(numeric(0))
            n <- max(1, round(5 * end / diff(lim)))
            b <- scales::breaks_extended(n + 1)(c(0, end))
            c(b[b > 0 & b < end * 0.8], end)
        }
        colour_scale <- scale_colour_gradientn(
            colours = pal,
            limits = lim,
            rescaler = rescale_diverging(0),
            oob = scales::oob_squish,
            name = bold_title(quote(-log[10] ~ "padj")),
            breaks = function(l) c(-rev(side_breaks(-l[1])), side_breaks(l[2])),
            labels = function(b) {
                out <- vapply(abs(signif(b, 3)), format, "", trim = TRUE)
                at_end <- (b <= lim[1] * 0.999 & sat[1]) | (b >= lim[2] * 0.999 & sat[2])
                out[at_end] <- paste0("\u2265", out[at_end])
                out
            },
            guide = colourbar
        )
    } else if (color.by == "direction") {
        df$colour_value <- factor(ifelse(df$direction == "up", "Up", "Down"), levels = c("Up", "Down"))
        colour_scale <- scale_colour_manual(
            values = c(Up = up_down[["up"]], Down = up_down[["down"]]),
            name = NULL,
            guide = guide_legend(override.aes = list(size = 2, alpha = 1))
        )
    } else {
        df$colour_value <- df$n_shared
        colour_scale <- scale_colour_gradientn(
            colours = resolve_continuous_palette(cols, FALSE),
            limits = c(1, max(2, length(cmp))),
            breaks = function(l) unique(round(scales::breaks_extended(4)(l))),
            name = "DE in n\ncomparisons"
        )
    }

    ns_df <- df[df$direction == "ns", , drop = FALSE]
    de_df <- df[df$direction != "ns", , drop = FALSE]
    # The most significant points last, on top
    de_df <- de_df[base::order(de_df$neg_log10_padj), , drop = FALSE]
    layer_fun <- points_layer_fun(raster, nrow(df), 5000, raster.dpi)
    # Positions along the comparison axis, values along the other
    map <- if (flip) aes(x = .data$logfc, y = .data$pos) else aes(x = .data$pos, y = .data$logfc)
    p <- ggplot(de_df, map)
    # Cutoffs as dashed lines, as in plot_volcano() (0 when there is none)
    cut_at <- if (logfc.cutoff > 0) c(-1, 1) * logfc.cutoff else 0
    p <- p + if (flip) {
        geom_vline(xintercept = cut_at, linetype = "dashed", colour = "grey40", linewidth = lw)
    } else {
        geom_hline(yintercept = cut_at, linetype = "dashed", colour = "grey40", linewidth = lw)
    }
    if (nrow(ns_df)) {
        # Not significant: grey, under the DE genes, as in plot_volcano()
        p <- p + layer_fun(geom_point(
            data = ns_df,
            colour = up_down[["ns"]],
            size = pt.size,
            alpha = alpha,
            shape = 16,
            stroke = 0
        ))
    }
    p <- p +
        layer_fun(geom_point(
            aes(colour = .data$colour_value),
            size = pt.size,
            alpha = alpha,
            shape = 16,
            stroke = 0
        )) +
        colour_scale

    if (show.counts) {
        counts <- stats::aggregate(
            cbind(up = direction == "up", down = direction == "down") ~ comparison,
            data = df,
            FUN = sum,
            drop = FALSE
        )
        counts$pos <- unname(pos[as.character(counts$comparison)])
        if (!is.null(grp_of)) counts$group <- grp_of[match(as.character(counts$comparison), cmp)]
        text_size <- pt_to_size(sz$small)
        p <- p + if (flip) {
            list(
                geom_text(data = counts, aes(x = vlim[2], y = .data$pos, label = .data$up),
                          inherit.aes = FALSE, hjust = 1, size = text_size, colour = up_down[["up"]]),
                geom_text(data = counts, aes(x = vlim[1], y = .data$pos, label = .data$down),
                          inherit.aes = FALSE, hjust = 0, size = text_size, colour = up_down[["down"]])
            )
        } else {
            list(
                geom_text(data = counts, aes(x = .data$pos, y = vlim[2], label = .data$up),
                          inherit.aes = FALSE, vjust = 1, size = text_size, colour = up_down[["up"]]),
                geom_text(data = counts, aes(x = .data$pos, y = vlim[1], label = .data$down),
                          inherit.aes = FALSE, vjust = 0, size = text_size, colour = up_down[["down"]])
            )
        }
    }
    if (any(df$labelled)) {
        check_installed("ggrepel", reason = "to label genes")
        # Labels of up-regulated genes stay beyond the upper cutoff and those
        # of down-regulated genes beyond the lower one, clear of the counts
        hi <- vlim[2] - room * unit
        lo <- vlim[1] + room * unit
        for (d in c("up", "down", "ns")) {
            lab_df <- df[df$labelled & df$direction == d, , drop = FALSE]
            if (nrow(lab_df) == 0) next
            span <- switch(d, up = c(logfc.cutoff, hi), down = c(lo, -logfc.cutoff), ns = c(lo, hi))
            repel <- list(
                data = lab_df,
                mapping = aes(label = .data$feature),
                size = pt_to_size(sz$small) * 0.9,
                fontface = if (italic.genes) "bold.italic" else "bold",
                colour = "black",
                max.overlaps = Inf,
                min.segment.length = 0,
                segment.colour = "grey30",
                segment.size = lw,
                box.padding = 0.2,
                max.time = 3,
                max.iter = 1e5,
                seed = seed
            )
            repel[[if (flip) "xlim" else "ylim"]] <- span
            p <- p + do.call(ggrepel::geom_text_repel, repel)
        }
    }

    # Comparison strips: a bar of each comparison's colour along its column
    # (above it) or row (left of it), its name beyond
    strips <- data.frame(fill = unname(cmp_cols[cmp]))
    strips[[if (flip) "y" else "x"]] <- unname(pos[cmp])
    strips[[if (flip) "yend" else "xend"]] <- unname(pos[cmp])
    if (!is.null(grp_of)) strips$group <- grp_of
    p <- p + group_bar_layer(strips, side = if (flip) "left" else "top", width = group.bar.width, half = 0.47)

    # Groups of comparisons: a bar along each group's panel and its name,
    # opposite the comparison names by default
    groups_side_abs <- if (!flip) {
        if (groups.side == "opposite") "bottom" else "top"
    } else {
        if (groups.side == "opposite") "right" else "left"
    }
    if (!is.null(grp_of)) {
        grp_lv <- levels(droplevels(grp_of))
        grp_bar_cols <- group_bar_colours(groups.cols, grp_lv, "groups.cols")
        gb <- data.frame(group = factor(grp_lv, levels = grp_lv), fill = unname(grp_bar_cols[grp_lv]))
        gb[[if (flip) "y" else "x"]] <- vapply(grp_lv, function(g) min(pos[cmp[grp_of == g]]), numeric(1))
        gb[[if (flip) "yend" else "xend"]] <- vapply(grp_lv, function(g) max(pos[cmp[grp_of == g]]), numeric(1))
        gap <- grid::unit(GROUP_BAR_GAP, "pt")
        # On the side of the comparison names: beyond the names
        offset <- if (groups_side_abs %in% c("bottom", "right")) {
            gap
        } else {
            name_room <- GROUP_BAR_GAP + group.bar.width + GROUP_BAR_TEXT_GAP + 1
            grid::unit(name_room, "pt") +
                label_extent(unname(cmp_text), sz$small, angle = if (upright) 0 else 45,
                             across = if (flip) "width" else "height") +
                gap
        }
        p <- p + group_bar_layer(gb, side = groups_side_abs, width = group.bar.width, offset = offset)
    }

    # Comparison axis: one unit per comparison; each panel spans its own
    n_per <- if (is.null(grp_of)) length(cmp) else as.vector(table(droplevels(grp_of)))
    if (is.null(grp_of)) {
        blank <- data.frame(pos = range(pos) + c(-0.5, 0.5), logfc = 0)
    } else {
        blank <- do.call(rbind, lapply(levels(droplevels(grp_of)), function(g) {
            data.frame(pos = range(pos[cmp[grp_of == g]]) + c(-0.5, 0.5), logfc = 0,
                       group = factor(g, levels = levels(droplevels(grp_of))))
        }))
    }
    name_margin <- group_bar_margin(if (flip) "left" else "top", group.bar.width)
    cmp_axis <- function(fun, position) {
        fun(breaks = unname(pos[cmp]), labels = unname(cmp_text[cmp]), position = position,
            expand = expansion(0))
    }
    p <- p + geom_blank(data = blank, aes(x = if (flip) .data$logfc else .data$pos,
                                          y = if (flip) .data$pos else .data$logfc), inherit.aes = FALSE)
    if (!flip) {
        p <- p +
            cmp_axis(scale_x_continuous, "top") +
            # oob_keep: values are never removed (no capping)
            scale_y_continuous(limits = vlim, expand = expansion(0), oob = scales::oob_keep) +
            labs(x = NULL, y = quote(log[2] ~ "fold change"))
    } else {
        p <- p +
            cmp_axis(scale_y_continuous, "left") +
            scale_x_continuous(limits = vlim, expand = expansion(0), oob = scales::oob_keep) +
            labs(x = quote(log[2] ~ "fold change"), y = NULL)
    }
    # A long last name, turned 45 degrees, reaches right of the panel
    right_margin <- if (!upright) {
        max(0.2, nchar(cmp_text[cmp[length(cmp)]]) * 0.55 * sz$small / 72 * 0.71 - cell.width / 2)
    } else {
        0.2
    }
    p <- p +
        theme_bioToolkit("axis", bg = bg) +
        (if (!flip) {
            panel_size_theme(n_per * cell.width, panel.height, panel.spacing = 0.1)
        } else {
            panel_size_theme(panel.width, n_per * cell.height, panel.spacing = 0.1)
        }) +
        coord_cartesian(clip = "off") +
        theme(
            legend.position = legend.position,
            legend.title = element_text(size = sz$small),
            plot.margin = margin(0.2, right_margin, 0.2, 0.2, "in"),
            strip.placement = "outside",
            strip.switch.pad.grid = unit(0, "pt"),
            strip.clip = "off",
            strip.text.x.top = element_text(size = sz$small, face = "bold",
                                            margin = group_bar_margin("top", group.bar.width)),
            strip.text.x.bottom = element_text(size = sz$small, face = "bold",
                                               margin = group_bar_margin("bottom", group.bar.width)),
            strip.text.y.left = element_text(size = sz$small, face = "bold", angle = 0, hjust = 1,
                                             margin = group_bar_margin("left", group.bar.width)),
            strip.text.y.right = element_text(size = sz$small, face = "bold", angle = 0, hjust = 0,
                                              margin = group_bar_margin("right", group.bar.width))
        )
    p <- p + if (!flip) {
        theme(
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x.top = element_text(
                angle = if (upright) 0 else 45,
                hjust = if (upright) 0.5 else 0,
                vjust = 0,
                margin = name_margin
            )
        )
    } else {
        theme(
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y.left = element_text(hjust = 1, margin = name_margin)
        )
    }
    if (!is.null(grp_of)) {
        p <- p + if (!flip) {
            facet_grid(cols = vars(.data$group), scales = "free_x", space = "free_x",
                       switch = if (groups_side_abs == "bottom") "x")
        } else {
            facet_grid(rows = vars(.data$group), scales = "free_y", space = "free_y",
                       switch = if (groups_side_abs == "left") "y")
        }
    }
    df$x <- df$pos
    src <- df[, intersect(c("group", "comparison", "feature", "x", "logfc", "padj", "direction",
                            "n_shared", "labelled"), names(df))]
    attr(p, "source_data") <- src
    p
}


#' Horizontal offsets of points within their column
#'
#' @param df Genes with `comparison`, `direction`, `logfc` and
#'           `neg_log10_padj`
#'
#' @param layout "sina", "pvalue" or "jitter"
#'
#' @param half Half the width of a column, in column units
#'
#' @return Offsets in [-half, half], along the rows of `df`
#'
#' @keywords internal
multi_volcano_offsets <- function(df,
                                  layout,
                                  half) {
    off <- numeric(nrow(df))
    key <- paste(df$comparison, df$direction)
    for (k in unique(key)) {
        i <- which(key == k)
        y <- df$logfc[i]
        off[i] <- switch(
            layout,
            jitter = stats::runif(length(i), -half, half),
            pvalue = {
                s <- df$neg_log10_padj[i]
                top <- max(s)
                if (top > 0) -half + 2 * half * s / top else rep(0, length(i))
            },
            sina = {
                # Width follows the density of genes at each fold change, the
                # widest point of each direction filling the column
                if (length(i) < 3 || stats::sd(y) == 0) {
                    stats::runif(length(i), -half, half) * 0.3
                } else {
                    d <- stats::density(y, n = 256)
                    w <- stats::approx(d$x, d$y, xout = y)$y
                    stats::runif(length(i), -1, 1) * half * w / max(d$y)
                }
            }
        )
    }
    off
}
