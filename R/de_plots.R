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
#   (diverging, symmetric), size = -log10 adjusted p-value. Significance
#   needs both cutoffs (padj and |log2FC|), so it is not readable from the
#   dot size alone: significant results are drawn opaque with a dark ring,
#   the others faded (default, "fade"), keeping their fold-change colour so
#   that no information is lost. One titled legend ("Significance") spells
#   out both cutoffs. "grey" and "hide" remain available. Replaces
#   `de_dot_plot()`.
# * `plot_multi_volcano()` (G4) shows *all* DE genes: one column of genes per
#   comparison, log2 fold change on the y-axis, a coloured band at 0 spanning
#   exactly +-`logfc.cutoff` (the non-significant zone). No capping of extreme
#   values (decided 2026-10-13).


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
#' @param logfc.limits Limits of the colour scale. If NULL (default),
#'                     symmetric around 0 at the largest absolute log2 fold
#'                     change of the significant dots. Values beyond get the
#'                     colour of the end
#'
#' @param cluster.features,cluster.comparisons Whether to order genes (within
#'        their groups) / comparisons (within their groups) by hierarchical
#'        clustering of their log2 fold changes, default is FALSE
#'
#' @param cols Colours of the diverging log2 fold change scale: NULL
#'             (default, "blue-magenta"), a palette name, or colours from low
#'             to high, e.g. `c("#2166AC", "white", "#B2182B")`
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
#' @details Genes not tested in a comparison have no dot. The source data
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
    flip = FALSE,
    italic.genes = TRUE,
    cell.width = 0.18,
    cell.height = 0.18,
    legend.position = c("right", "left", "top", "bottom", "none"),
    bg = "transparent"
) {
    method <- match.arg(method)
    significance <- match.arg(significance)
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
        v <- abs(df$logfc[drawn & df$significant])
        if (length(v) == 0) v <- abs(df$logfc[drawn])
        m <- max(c(v, 0), na.rm = TRUE)
        if (m == 0) m <- 1
        logfc.limits <- c(-m, m)
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
    feat_labels <- stats::setNames(ft$feature, ft$id)
    cmp_labels <- stats::setNames(cmp, cmp)
    italic <- stats::setNames(rep(italic.genes, length(feat_labels)), names(feat_labels))

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
        cell.width = cell.width,
        cell.height = cell.height,
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
#'              column of DE genes per comparison, at their log2 fold change,
#'              with a coloured band at 0 spanning exactly
#'              +-`logfc.cutoff` (the non-significant zone), the numbers of up-
#'              and down-regulated genes above and below each column, and the
#'              top genes labelled. An improved version of the plot made
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
#' @param show.ns Whether to draw non-significant genes too, faint and under
#'                the band, default is FALSE
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
#' @param band.cols Colours of the bands, one per comparison: NULL (default;
#'                  the default group palette), a qualitative palette name, or
#'                  colours named by comparison
#'
#' @param band.labels Display names of the comparisons, named by comparison.
#'                    Names are written in the bands when they fit, otherwise
#'                    under the columns
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
#' @param raster Whether to rasterize the points (labels, bands and text stay
#'               vector): NULL (default) rasterizes above 5,000 points, or
#'               TRUE / FALSE
#'
#' @param raster.dpi Resolution of rasterized points, default is 1000
#'
#' @param cell.width Width of one column in inches, default is 0.6
#'
#' @param panel.height Height of the panel in inches, default is 3.5
#'
#' @param legend.position Position of the legend, default is "right"
#'
#' @param bg Background colour, default is "transparent"
#'
#' @param seed Random seed for the spread of points and label placement,
#'             default is 42
#'
#' @details The source data (\code{\link{export_source_data}}) hold one row
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
    show.ns = FALSE,
    layout = c("sina", "pvalue", "jitter"),
    color.by = c("padj", "direction", "sharing"),
    cols = NULL,
    band.cols = NULL,
    band.labels = NULL,
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
    panel.height = 3.5,
    legend.position = "right",
    bg = "transparent",
    seed = 42
) {
    method <- match.arg(method)
    layout <- match.arg(layout)
    color.by <- match.arg(color.by)
    label <- match.arg(label)
    label.rank <- match.arg(label.rank)
    order <- match.arg(order)
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

    df <- res[is_de | show.ns, , drop = FALSE]
    if (!any(df$direction != "ns")) warning("No DE genes at these cutoffs", call. = FALSE)
    df$comparison <- factor(as.character(df$comparison), levels = cmp)
    # Column position: within the group's panel when grouped
    if (is.null(grp_of)) {
        pos <- stats::setNames(seq_along(cmp), cmp)
    } else {
        pos <- stats::setNames(stats::ave(seq_along(cmp), grp_of, FUN = seq_along), cmp)
    }
    half <- 0.4
    df$x <- pos[as.character(df$comparison)] +
        with_seed(seed, multi_volcano_offsets(df, layout, half))

    # y limits: all values (no capping), with room in inches for the counts
    y_rng <- range(c(df$logfc, -logfc.cutoff, logfc.cutoff, 0))
    room <- if (show.counts) 0.22 else 0.05
    unit <- diff(y_rng) / (panel.height - 2 * room)
    ylim <- y_rng + c(-1, 1) * room * unit
    # Band: +-logfc.cutoff, at least tall enough for its label
    band_text <- if (is.null(band.labels)) stats::setNames(cmp, cmp) else {
        lab <- stats::setNames(cmp, cmp)
        lab[intersect(names(band.labels), cmp)] <- band.labels[intersect(names(band.labels), cmp)]
        lab
    }
    band_h <- max(logfc.cutoff, 0.6 * sz$small / 72 * unit)
    fits <- all(nchar(band_text) * 0.55 * sz$small / 72 <= 0.92 * cell.width) &&
        2 * band_h / unit >= 1.1 * sz$small / 72
    band_cols <- group_colours(band.cols, factor(cmp, levels = cmp), arg = "band.cols")$cols
    band <- data.frame(comparison = factor(cmp, levels = cmp), x = unname(pos[cmp]),
                       fill = unname(band_cols[cmp]), label = unname(band_text[cmp]))
    if (!is.null(grp_of)) band$group <- grp_of
    lum <- colSums(grDevices::col2rgb(band$fill) * c(0.299, 0.587, 0.114))
    band$text_col <- ifelse(lum < 140, "white", "black")

    # Labels: top genes per column and direction, plus selected genes
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
    if (color.by == "padj") {
        df$colour_value <- ifelse(df$direction == "up", 1, -1) * df$neg_log10_padj
        # The colour saturates at the 95th percentile of the DE genes' -log10
        # padj: a few p-values near 0 (common with per-cell tests) would
        # otherwise leave every other gene pale. Positions are not capped
        v <- abs(df$colour_value[df$direction != "ns"])
        m <- max(if (length(v)) stats::quantile(v, 0.95, names = FALSE) else 0, -log10(padj.cutoff) * 2)
        saturated <- any(v > m)
        pal <- if (is.null(cols)) bt_palette("blue-magenta")[c(1:4, 8:11)] else resolve_continuous_palette(cols, TRUE)
        colour_scale <- scale_colour_gradientn(
            colours = pal, limits = c(-m, m), oob = scales::oob_squish, name = quote(-log[10] ~ "padj"),
            labels = function(b) {
                out <- format(abs(signif(b, 3)), trim = TRUE)
                if (saturated) out[abs(b) >= m * 0.999] <- paste0("\u2265", out[abs(b) >= m * 0.999])
                out
            },
            breaks = function(l) {
                b <- scales::breaks_extended(5)(c(0, l[2]))
                b <- b[b > 0 & b < l[2] * 0.9]
                c(-rev(c(b, l[2])), c(b, l[2]))
            },
            guide = guide_colourbar(theme = theme(
                legend.key.height = unit(1.2, "in"), legend.key.width = unit(0.12, "in"),
                legend.frame = element_rect(colour = "black", linewidth = 0.3),
                legend.ticks = element_line(colour = "black", linewidth = 0.3)
            ))
        )
    } else if (color.by == "direction") {
        df$colour_value <- factor(ifelse(df$direction == "up", "Up", "Down"), levels = c("Up", "Down"))
        colour_scale <- scale_colour_manual(values = c(Up = up_down[["up"]], Down = up_down[["down"]]), name = NULL,
                                            guide = guide_legend(override.aes = list(size = 2, alpha = 1)))
    } else {
        df$colour_value <- df$n_shared
        colour_scale <- scale_colour_gradientn(
            colours = resolve_continuous_palette(cols, FALSE), limits = c(1, max(2, length(cmp))),
            breaks = function(l) unique(round(scales::breaks_extended(4)(l))),
            name = "DE in n\ncomparisons"
        )
    }

    ns_df <- df[df$direction == "ns", , drop = FALSE]
    de_df <- df[df$direction != "ns", , drop = FALSE]
    # The most significant points last, on top
    de_df <- de_df[base::order(de_df$neg_log10_padj), , drop = FALSE]
    layer_fun <- points_layer_fun(raster, nrow(df), 5000, raster.dpi)
    p <- ggplot(de_df, aes(x = .data$x, y = .data$logfc))
    if (nrow(ns_df)) {
        p <- p + layer_fun(geom_point(data = ns_df, colour = "grey80", size = pt.size * 0.8, alpha = 0.5,
                                      shape = 16, stroke = 0))
    }
    p <- p +
        geom_rect(data = band, aes(xmin = .data$x - 0.47, xmax = .data$x + 0.47, ymin = -band_h, ymax = band_h),
                  fill = band$fill, inherit.aes = FALSE, colour = NA) +
        layer_fun(geom_point(aes(colour = .data$colour_value), size = pt.size, alpha = alpha,
                             shape = 16, stroke = 0)) +
        colour_scale
    if (fits) {
        p <- p + geom_text(data = band, aes(x = .data$x, y = 0, label = .data$label), inherit.aes = FALSE,
                           colour = band$text_col, size = pt_to_size(sz$small))
    }
    if (show.counts) {
        counts <- stats::aggregate(cbind(up = direction == "up", down = direction == "down") ~ comparison,
                                   data = df, FUN = sum, drop = FALSE)
        counts$x <- unname(pos[as.character(counts$comparison)])
        if (!is.null(grp_of)) counts$group <- grp_of[match(as.character(counts$comparison), cmp)]
        p <- p +
            geom_text(data = counts, aes(x = .data$x, y = ylim[2], label = .data$up), inherit.aes = FALSE,
                      vjust = 1, size = pt_to_size(sz$small), colour = up_down[["up"]]) +
            geom_text(data = counts, aes(x = .data$x, y = ylim[1], label = .data$down), inherit.aes = FALSE,
                      vjust = 0, size = pt_to_size(sz$small), colour = up_down[["down"]])
    }
    if (any(df$labelled)) {
        check_installed("ggrepel", reason = "to label genes")
        # Labels of up-regulated genes stay above the band and those of
        # down-regulated genes below it, clear of the counts
        top <- ylim[2] - room * unit
        bottom <- ylim[1] + room * unit
        for (d in c("up", "down", "ns")) {
            lab_df <- df[df$labelled & df$direction == d, , drop = FALSE]
            if (nrow(lab_df) == 0) next
            p <- p + ggrepel::geom_text_repel(
                data = lab_df, aes(label = .data$feature),
                size = pt_to_size(sz$small) * 0.9,
                fontface = if (italic.genes) "bold.italic" else "bold",
                colour = "black", max.overlaps = Inf, min.segment.length = 0,
                segment.colour = "grey30", segment.size = lw, box.padding = 0.2,
                max.time = 3, max.iter = 1e5, seed = seed,
                ylim = switch(d, up = c(band_h, top), down = c(bottom, -band_h), ns = c(bottom, top))
            )
        }
    }

    # x: one unit per column; the panel is as wide as its columns
    n_cols <- if (is.null(grp_of)) length(cmp) else as.vector(table(droplevels(grp_of)))
    blank <- data.frame(x = c(0.5, if (is.null(grp_of)) length(cmp) + 0.5 else n_cols + 0.5), logfc = 0)
    if (!is.null(grp_of)) {
        blank <- data.frame(x = c(rep(0.5, length(n_cols)), n_cols + 0.5), logfc = 0,
                            group = factor(rep(levels(droplevels(grp_of)), 2), levels = levels(droplevels(grp_of))))
    }
    p <- p +
        geom_blank(data = blank) +
        # Column names are in the bands when they fit; otherwise under the
        # columns (set below)
        scale_x_continuous(
            breaks = if (fits || !is.null(grp_of)) NULL else seq_along(cmp),
            labels = if (fits || !is.null(grp_of)) NULL else unname(band_text[cmp]),
            expand = expansion(0)
        ) +
        scale_y_continuous(limits = ylim, expand = expansion(0)) +
        labs(x = NULL, y = quote(log[2] ~ "fold change")) +
        theme_bioToolkit("axis", bg = bg) +
        panel_size_theme(n_cols * cell.width, panel.height, panel.spacing = 0.1) +
        theme(
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            # A long first column name, rotated, reaches left of the y-axis
            plot.margin = margin(0.2, 0.2, 0.2,
                                 max(0.2, nchar(band_text[cmp[1]]) * 0.55 * sz$small / 72 * 0.71 -
                                         cell.width / 2 - 0.45), "in"),
            legend.position = legend.position,
            legend.title = element_text(size = sz$small)
        )
    if (!fits && !is.null(grp_of)) {
        # Free x scales restart at 1 in each group's panel, so the names are
        # drawn as text under each panel instead of axis labels
        labs_df <- data.frame(x = unname(pos[cmp]), group = grp_of, label = unname(band_text[cmp]))
        p <- p + geom_text(data = labs_df, aes(x = .data$x, y = ylim[1], label = .data$label),
                           inherit.aes = FALSE, angle = 45, hjust = 1, vjust = 1,
                           nudge_y = -0.03 * unit, size = pt_to_size(sz$small)) +
            coord_cartesian(clip = "off") +
            theme(plot.margin = margin(0.2, 0.2, max(nchar(band_text)) * 0.45 * sz$small / 72 + 0.2, 0.2, "in"))
    }
    if (!is.null(grp_of)) {
        p <- p + facet_grid(cols = vars(.data$group), scales = "free_x", space = "free_x") +
            theme(strip.text = element_text(size = sz$small, face = "bold"))
    }
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
