# logFC-logFC quadrant plot of two DE comparisons (roadmap item B5)
#
# Design
# ------
# `plot_fc_quadrant()` compares two DE results feature by feature: the log2
# fold change in one comparison on x, in the other on y. Typical uses: the
# same contrast at two time points or in two tissues, a knockout vs a drug,
# or RNA vs chromatin (ATAC / ChIP) changes of the same genes.
#
# * Both tables are read with `standardize_de()` (F5), each with its own
#   detected method, so e.g. an edgeR RNA table can be compared with a DESeq2
#   ATAC table. Tables already standardized are used as they are.
# * Features are joined by name (inner join, reported in a message). For
#   peaks (ATAC / ChIP), `map.y` maps the features of `de.y` to genes; many
#   peaks can map to one gene, and `map.summary` picks one peak per gene: the
#   most significant one ("min_padj", default), whose fold change is the most
#   reliable, or the one with the largest absolute fold change ("max_abs"),
#   which can be a noisy, weakly supported peak. Averaging fold changes was
#   not offered: peaks of one gene can change in opposite directions, and
#   their mean has no p-value to decide significance.
# * Five classes, drawn in this order so that the informative ones are on
#   top: neither significant (grey), significant in x only, in y only (light
#   colours), in both with the same sign (concordant) or opposite signs
#   (discordant), slightly larger and in deep colours. The defaults follow
#   one idea: x only is blue, y only is red, both in the same direction is
#   purple (the two mixed), opposite directions is green (from Paul Tol's
#   colour-blind safe schemes).
# * Counts of the features significant in both in each quadrant, in the
#   corners and in the colour of their class; the correlation of the fold
#   changes (Spearman by default: robust to the long tails of fold changes)
#   and its n as a subtitle, outside the panel, where it never covers points.
# * The same symmetric range on both axes by default (`equal.axes`), so that
#   the y = x diagonal (dashed) is the line of equal change and a slope
#   flatter than it reads as a weaker effect in y. Square panel (F2). No
#   capping: every feature is drawn at its value.
# * Labels: the `label.n` features farthest from the origin in each quadrant,
#   among those significant in both (so discordant genes are named too), or
#   only `label.features`; bold italic as in the volcano.


#' @title Quadrant plot of the log fold changes of two DE comparisons
#'
#' @description Plots, for every feature measured in both comparisons, its
#'              log2 fold change in one (x) against the other (y), coloured by
#'              where it is significant: in both, with the same or opposite
#'              signs, in one only, or in neither. Shows how far two
#'              comparisons agree (e.g. two time points, two tissues, or RNA
#'              and chromatin), with counts of features significant in both
#'              per quadrant, the correlation of the fold changes, the y = x
#'              diagonal and labels of the top features. Square panel of a
#'              fixed size; save it with \code{\link{save_plot}}; the joined
#'              table behind it is returned by \code{\link{export_source_data}}.
#'
#' @param de.x,de.y DE results of the comparisons on x and y: data frames
#'                  (or DESeq2 `DESeqResults` / edgeR `TopTags`) read with
#'                  \code{\link{standardize_de}}, each with its own detected
#'                  method; or tables already returned by it. For tables whose
#'                  columns must be named and differ between the two,
#'                  standardize each first
#'
#' @param name.x,name.y Names of the comparisons, used in the axis titles and
#'                      the legend. Default "x" and "y"
#'
#' @param method,feature.col,logfc.col,p.col,padj.col How to read both
#'        tables, as in \code{\link{standardize_de}}. Default: detected
#'
#' @param map.y Mapping of the features of `de.y` to the features of `de.x`,
#'              e.g. peaks to genes (from a peak annotation such as
#'              ChIPseeker's, filtered as wanted, e.g. to promoter peaks): a
#'              data frame whose first column holds features of `de.y` and
#'              second the genes, or a character vector of genes named by
#'              feature. Features of `de.y` without a gene are left out.
#'              Default is NULL (features are matched by name)
#'
#' @param map.summary With `map.y`, how one feature (peak) is chosen per gene
#'                    when several map to it: "min_padj" (default; the
#'                    smallest adjusted p-value, ties broken by the largest
#'                    absolute fold change) or "max_abs" (the largest absolute
#'                    fold change)
#'
#' @param padj.cutoff Adjusted p-value below which a feature is significant,
#'                    default 0.05. One value for both comparisons, or two
#'                    (x, y)
#'
#' @param logfc.cutoff Absolute log2 fold change a significant feature must
#'                     exceed, default 0. One value, or two (x, y)
#'
#' @param categories Which classes are coloured: any of "both" (significant
#'                   in both, concordant and discordant), "x_only" and
#'                   "y_only". The others are drawn in grey with the features
#'                   significant in neither. Default all three
#'
#' @param show.counts Whether to print the number of features significant in
#'                    both in each quadrant, in its corner. Default TRUE
#'
#' @param correlation Correlation of the fold changes shown above the panel:
#'                    "spearman" (default), "pearson" or "none"
#'
#' @param cor.on Features the correlation is computed on: "all" (default;
#'               every feature drawn) or "significant" (significant in x or
#'               y)
#'
#' @param diagonal Whether to draw the y = x line (dashed), default TRUE
#'
#' @param equal.axes Whether x and y have the same range, symmetric around 0
#'                   (default TRUE), so that the diagonal is the line of equal
#'                   change. If FALSE, each axis is symmetric around 0 on its
#'                   own range
#'
#' @param label Which features to label: "top" (default; the `label.n`
#'              features farthest from the origin in each quadrant, among
#'              those significant in both, plus any `label.features`),
#'              "selected" (only `label.features`) or "none"
#'
#' @param label.n Number of features labelled per quadrant with
#'                `label = "top"`, default 5
#'
#' @param label.features Features to label, matched exactly (the case of
#'                       gene symbols matters)
#'
#' @param label.size Size of the labels in mm. If NULL (default), the small
#'                   text size of the theme (9 pt)
#'
#' @param italic.genes Whether labels are gene names, set in italics. Default
#'                     TRUE. Labels are bold either way
#'
#' @param cols Colours of the classes: a vector named by any of
#'             `concordant`, `discordant`, `x_only`, `y_only` and `ns` (the
#'             others keep their default), or five colours in that order.
#'             NULL (default) is `c(concordant = "#762A83", discordant =
#'             "#1B7837", x_only = "#4393C3", y_only = "#D6604D", ns =
#'             "grey80")`
#'
#' @param pt.size Size of the points, default 1; features significant in both
#'                are drawn 1.4 times as large
#'
#' @param alpha Opacity of the points, default 0.7
#'
#' @param raster Whether to rasterize the points while keeping text and lines
#'               as vectors: NULL (default) rasterizes when more than 5,000
#'               features are drawn, or TRUE / FALSE
#'
#' @param raster.dpi Resolution of the rasterized points, default 1000
#'
#' @param panel.width Width of the panel in inches, default 3
#'
#' @param panel.height Height of the panel in inches. Default is
#'                     `panel.width` (a square panel)
#'
#' @param legend.position Position of the legend, default "right"
#'
#' @param bg Background colour, default is "transparent"
#'
#' @param seed Random seed for the placement of labels, default is 42
#'
#' @details
#' A feature is significant in a comparison when its adjusted p-value is
#' below `padj.cutoff` and its absolute log2 fold change above
#' `logfc.cutoff`; features with an NA adjusted p-value (e.g. DESeq2
#' independent filtering) are not significant. Features without a fold change
#' in either table are not drawn.
#'
#' The correlation describes the agreement of all fold changes, most of them
#' small and noisy; `cor.on = "significant"` restricts it to the features
#' changed in at least one comparison. Neither is a test of whether the two
#' comparisons share biology.
#'
#' The source data (\code{\link{export_source_data}}) hold one row per
#' drawn feature: `feature`, `logfc_x`, `padj_x`, `logfc_y`, `padj_y`,
#' `category`, `labelled` and, with `map.y`, the feature of `de.y` chosen
#' (`feature_y`) and the number mapped to the gene (`n_mapped`).
#'
#' @return A ggplot object
#'
#' @examples
#' set.seed(1)
#' genes <- paste0("Gene", 1:3000)
#' shared <- rnorm(3000, sd = 1.2)
#' de1 <- data.frame(logFC = shared + rnorm(3000, sd = 0.5), logCPM = 5, row.names = genes)
#' de2 <- data.frame(logFC = 0.7 * shared + rnorm(3000, sd = 0.5), logCPM = 5, row.names = genes)
#' de1$PValue <- 2 * pnorm(-abs(de1$logFC) / 0.45)
#' de2$PValue <- 2 * pnorm(-abs(de2$logFC) / 0.45)
#' de1$FDR <- p.adjust(de1$PValue, "BH")
#' de2$FDR <- p.adjust(de2$PValue, "BH")
#' plot_fc_quadrant(de1, de2, name.x = "1 week", name.y = "6 weeks")
#' \dontrun{
#' # RNA vs ATAC: peaks mapped to genes, most significant peak per gene
#' peak2gene <- data.frame(peak = anno$peak_id, gene = anno$SYMBOL)
#' p <- plot_fc_quadrant(rna_de, atac_de, name.x = "RNA", name.y = "ATAC", map.y = peak2gene)
#' save_plot(p, "quadrant.pdf")
#' }
#'
#' @seealso \code{\link{plot_volcano}}, \code{\link{standardize_de}},
#'          \code{\link{save_plot}}
#'
#' @import ggplot2
#' @export
#'
plot_fc_quadrant <- function(
    de.x,
    de.y,
    name.x = "x",
    name.y = "y",
    method = c("auto", "edgeR", "DESeq2", "limma", "seurat"),
    feature.col = NULL,
    logfc.col = NULL,
    p.col = NULL,
    padj.col = NULL,
    map.y = NULL,
    map.summary = c("min_padj", "max_abs"),
    padj.cutoff = 0.05,
    logfc.cutoff = 0,
    categories = c("both", "x_only", "y_only"),
    show.counts = TRUE,
    correlation = c("spearman", "pearson", "none"),
    cor.on = c("all", "significant"),
    diagonal = TRUE,
    equal.axes = TRUE,
    label = c("top", "selected", "none"),
    label.n = 5,
    label.features = NULL,
    label.size = NULL,
    italic.genes = TRUE,
    cols = NULL,
    pt.size = 1,
    alpha = 0.7,
    raster = NULL,
    raster.dpi = 1000,
    panel.width = 3,
    panel.height = panel.width,
    legend.position = "right",
    bg = "transparent",
    seed = 42
) {
    method <- match.arg(method)
    map.summary <- match.arg(map.summary)
    correlation <- match.arg(correlation)
    cor.on <- match.arg(cor.on)
    label <- match.arg(label)
    categories <- match.arg(categories, several.ok = TRUE)
    padj.cutoff <- rep_len(padj.cutoff, 2)
    logfc.cutoff <- rep_len(logfc.cutoff, 2)
    for (i in 1:2) check_de_cutoffs(padj.cutoff[i], logfc.cutoff[i])
    cols <- quadrant_colours(cols)

    read <- function(d, what) {
        # A table already standardized is used as it is
        if (is.data.frame(d) && all(c("feature", "logfc", "padj") %in% names(d)) && !is.null(attr(d, "method"))) {
            out <- d
        } else {
            out <- standardize_de(d, method = method, feature.col = feature.col, logfc.col = logfc.col,
                                  p.col = p.col, padj.col = padj.col)
        }
        if (!is.null(out$comparison)) stop("`", what, "` must hold one comparison", call. = FALSE)
        out[c("feature", "logfc", "padj")]
    }
    dx <- read(de.x, "de.x")
    dy <- read(de.y, "de.y")
    if (!is.null(map.y)) dy <- map_features(dy, map.y, map.summary)
    for (nm in c("dx", "dy")) {
        d <- get(nm)
        if (anyDuplicated(d$feature)) {
            stop("Features appear more than once in `", if (nm == "dx") "de.x" else "de.y", "`",
                 if (nm == "dy" && is.null(map.y)) "; to collapse many peaks per gene, use `map.y`",
                 call. = FALSE)
        }
    }

    # Features measured in both, with a fold change in both
    df <- merge(dx, dy, by = "feature", suffixes = c("_x", "_y"))
    df <- df[!is.na(df$logfc_x) & !is.na(df$logfc_y), , drop = FALSE]
    message(nrow(df), " feature(s) in both comparisons (", nrow(dx), " in x, ", nrow(dy), " in y)")
    if (nrow(df) < 3) stop("Fewer than 3 features are shared by the two comparisons", call. = FALSE)

    dir_x <- de_direction(df$logfc_x, df$padj_x, padj.cutoff[1], logfc.cutoff[1])
    dir_y <- de_direction(df$logfc_y, df$padj_y, padj.cutoff[2], logfc.cutoff[2])
    sig_x <- dir_x != "ns"
    sig_y <- dir_y != "ns"
    cls <- ifelse(sig_x & sig_y, ifelse(dir_x == dir_y, "concordant", "discordant"),
                  ifelse(sig_x, "x_only", ifelse(sig_y, "y_only", "ns")))
    # Classes not asked for are drawn as not significant
    shown <- c(if ("both" %in% categories) c("concordant", "discordant"), intersect(c("x_only", "y_only"), categories))
    cls[!(cls %in% shown)] <- "ns"
    class_labels <- c(
        concordant = "Both, same direction",
        discordant = "Both, opposite",
        x_only = paste(name.x, "only"),
        y_only = paste(name.y, "only"),
        ns = "Neither"
    )
    lv <- c(shown, "ns")
    df$category <- factor(class_labels[cls], levels = class_labels[lv])
    df$labelled <- quadrant_label_flags(df, cls, label, label.n, label.features)
    # Drawn from the least to the most informative class
    draw_order <- c("ns", "x_only", "y_only", "discordant", "concordant")
    # Significant in x or y, whether or not that class is coloured (for the
    # correlation on significant features)
    df$.any_sig <- sig_x | sig_y
    ord <- order(match(cls, draw_order))
    df <- df[ord, , drop = FALSE]
    cls <- cls[ord]
    rownames(df) <- NULL
    df$.size <- ifelse(cls %in% c("concordant", "discordant"), pt.size * 1.4, pt.size)

    sz <- theme_sizes()
    lw <- pt_to_linewidth(sz$line)
    if (is.null(label.size)) label.size <- pt_to_size(sz$small)
    m_x <- max(abs(df$logfc_x))
    m_y <- max(abs(df$logfc_y))
    if (equal.axes) m_x <- m_y <- max(m_x, m_y)
    xlim <- c(-1, 1) * m_x
    ylim <- c(-1, 1) * m_y

    layer_fun <- points_layer_fun(raster, nrow(df), 5000, raster.dpi)
    p <- ggplot(df, aes(x = .data$logfc_x, y = .data$logfc_y)) +
        geom_hline(yintercept = 0, colour = "grey70", linewidth = lw) +
        geom_vline(xintercept = 0, colour = "grey70", linewidth = lw)
    if (diagonal) {
        # Stops short of the corners (at 88% of the shorter half-axis), where
        # the counts are printed
        d_end <- 0.88 * min(m_x, m_y) * 1.06
        p <- p + annotate("segment", x = -d_end, y = -d_end, xend = d_end, yend = d_end, linetype = "dashed",
                          colour = "grey45", linewidth = lw)
    }
    p <- p +
        layer_fun(geom_point(aes(colour = .data$category, size = .data$.size), alpha = alpha, shape = 16,
                             stroke = 0)) +
        scale_size_identity() +
        scale_colour_manual(values = stats::setNames(cols[lv], class_labels[lv]), name = NULL, drop = TRUE) +
        scale_x_continuous(limits = xlim, expand = expansion(mult = 0.06)) +
        scale_y_continuous(limits = ylim, expand = expansion(mult = 0.06)) +
        labs(x = bquote(log[2] ~ "fold change," ~ .(name.x)), y = bquote(log[2] ~ "fold change," ~ .(name.y))) +
        theme_bioToolkit("axis", bg = bg) +
        panel_size_theme(panel.width, panel.height) +
        theme(legend.position = legend.position,
              plot.subtitle = element_text(size = sz$small, colour = "grey25", hjust = 0,
                                           margin = margin(b = 4))) +
        guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1)))

    if (show.counts && "both" %in% categories) {
        both <- cls %in% c("concordant", "discordant")
        q <- quadrant_counts(df$logfc_x[both], df$logfc_y[both])
        count_size <- pt_to_size(sz$base)
        # Same-sign quadrants (top right, bottom left) in the concordant
        # colour, opposite-sign ones in the discordant colour
        p <- p + annotate(
            "text",
            x = c(Inf, -Inf, -Inf, Inf),
            y = c(Inf, Inf, -Inf, -Inf),
            label = format(q, big.mark = ","),
            hjust = c(1.3, -0.3, -0.3, 1.3),
            vjust = c(1.5, 1.5, -0.6, -0.6),
            size = count_size,
            fontface = "bold",
            colour = cols[c("concordant", "discordant", "concordant", "discordant")]
        )
    }

    if (correlation != "none") {
        use <- if (cor.on == "all") rep(TRUE, nrow(df)) else df$.any_sig
        n_use <- sum(use)
        if (n_use >= 3) {
            r <- stats::cor(df$logfc_x[use], df$logfc_y[use], method = correlation)
            r_txt <- formatC(r, format = "f", digits = 2)
            n_txt <- format(n_use, big.mark = ",")
            sub <- if (correlation == "spearman") {
                bquote("Spearman" ~ rho == .(r_txt) * "," ~ italic(n) == .(n_txt))
            } else {
                bquote("Pearson" ~ italic(r) == .(r_txt) * "," ~ italic(n) == .(n_txt))
            }
            p <- p + labs(subtitle = sub)
        } else {
            message("Fewer than 3 features for the correlation; not shown")
        }
    }

    if (any(df$labelled)) {
        check_installed("ggrepel", reason = "to label features")
        p <- p + ggrepel::geom_text_repel(
            data = df[df$labelled, , drop = FALSE],
            aes(label = .data$feature),
            size = label.size,
            fontface = if (italic.genes) "bold.italic" else "bold",
            colour = "black",
            max.overlaps = Inf,
            min.segment.length = 0,
            segment.colour = "grey30",
            segment.size = lw,
            box.padding = 0.3,
            max.time = 3,
            max.iter = 1e5,
            seed = seed
        )
    }
    src <- df[setdiff(names(df), c(".size", ".any_sig"))]
    attr(p, "source_data") <- src
    p
}


#' Colours of the classes of a quadrant plot
#'
#' @param cols NULL, colours named by class (the others keep their default),
#'             or five unnamed colours in the order concordant, discordant,
#'             x_only, y_only, ns
#'
#' @return Five colours named by class
#'
#' @keywords internal
quadrant_colours <- function(cols) {
    # Paul Tol's PRGn ends (deep purple, deep green) for the features
    # significant in both, and the blue and red of his BuRd for one only:
    # blue (x) and red (y) mix into purple (both, same direction)
    out <- c(concordant = "#762A83", discordant = "#1B7837", x_only = "#4393C3", y_only = "#D6604D",
             ns = "grey80")
    if (is.null(cols)) return(out)
    if (!all(is_colour(cols))) stop("`cols` must be colours", call. = FALSE)
    if (is.null(names(cols))) {
        if (length(cols) != 5) {
            stop("Unnamed `cols` must be five colours: concordant, discordant, x_only, y_only, ns", call. = FALSE)
        }
        names(cols) <- names(out)
    }
    unknown <- setdiff(names(cols), names(out))
    if (length(unknown)) {
        stop("`cols` must be named by: ", paste(names(out), collapse = ", "), "; not: ",
             paste(unknown, collapse = ", "), call. = FALSE)
    }
    out[names(cols)] <- cols
    out
}


#' Map features of a DE table (e.g. peaks) to genes, one per gene
#'
#' @param d Standardized DE table with `feature`, `logfc`, `padj`
#'
#' @param map A data frame (first column: features, second: genes) or genes
#'            named by feature
#'
#' @param summary "min_padj" or "max_abs"
#'
#' @return `d` with one row per gene: `feature` is the gene, `feature_y` the
#'         feature chosen and `n_mapped` the number of features of the gene
#'
#' @keywords internal
map_features <- function(d,
                         map,
                         summary) {
    if (is.data.frame(map)) {
        if (ncol(map) < 2) stop("`map.y` must have two columns: features, genes", call. = FALSE)
        map <- stats::setNames(as.character(map[[2]]), as.character(map[[1]]))
    }
    if (!is.character(map) || is.null(names(map))) {
        stop("`map.y` must be a data frame (features, genes) or genes named by feature", call. = FALSE)
    }
    map <- map[!is.na(map) & nzchar(map)]
    gene <- unname(map[d$feature])
    keep <- !is.na(gene)
    if (!any(keep)) stop("No feature of `de.y` is in `map.y`", call. = FALSE)
    d <- d[keep, , drop = FALSE]
    d$feature_y <- d$feature
    d$feature <- gene[keep]
    # Best feature first, then the first row of each gene is kept
    key <- if (summary == "min_padj") {
        order(d$feature, ifelse(is.na(d$padj), Inf, d$padj), -abs(d$logfc))
    } else {
        order(d$feature, -abs(d$logfc), ifelse(is.na(d$padj), Inf, d$padj))
    }
    d <- d[key, , drop = FALSE]
    n <- table(d$feature)
    d <- d[!duplicated(d$feature), , drop = FALSE]
    d$n_mapped <- as.integer(n[d$feature])
    message(sum(keep), " feature(s) of `de.y` mapped to ", nrow(d), " gene(s); one per gene by ",
            if (summary == "min_padj") "smallest adjusted p-value" else "largest absolute fold change")
    d
}


#' Which features of a quadrant plot are labelled
#'
#' @param df Joined table with `feature`, `logfc_x`, `logfc_y`
#'
#' @param cls Class of each row
#'
#' @param label "top", "selected" or "none"
#'
#' @param label.n Features per quadrant
#'
#' @param label.features Features to label, or NULL
#'
#' @return A logical vector along the rows of `df`
#'
#' @keywords internal
quadrant_label_flags <- function(df,
                                 cls,
                                 label,
                                 label.n,
                                 label.features) {
    flags <- rep(FALSE, nrow(df))
    if (label == "none") return(flags)
    if (label == "selected" && length(label.features) == 0) {
        stop("`label = \"selected\"` needs `label.features`", call. = FALSE)
    }
    if (label == "top") {
        if (!(is.numeric(label.n) && length(label.n) == 1 && label.n >= 0)) {
            stop("`label.n` must be a single number >= 0", call. = FALSE)
        }
        both <- cls %in% c("concordant", "discordant")
        quadrant <- paste(sign(df$logfc_x), sign(df$logfc_y))
        dist <- sqrt(df$logfc_x^2 + df$logfc_y^2)
        for (q in unique(quadrant[both])) {
            idx <- which(both & quadrant == q)
            flags[idx[utils::head(order(-dist[idx]), label.n)]] <- TRUE
        }
    }
    if (length(label.features)) {
        label.features <- as.character(label.features)
        missing <- setdiff(label.features, df$feature)
        if (length(missing)) {
            near <- df$feature[tolower(df$feature) %in% tolower(missing)]
            warning("`label.features` not found among the drawn features: ", paste(missing, collapse = ", "),
                    if (length(near)) paste0(" (matching only in case: ", paste(unique(near), collapse = ", "), ")"),
                    call. = FALSE)
        }
        flags <- flags | df$feature %in% label.features
    }
    flags
}


#' Number of points in each quadrant
#'
#' Points on an axis (a fold change of exactly 0) belong to no quadrant.
#'
#' @param x,y Coordinates
#'
#' @return Counts in the order top right, top left, bottom left, bottom
#'         right
#'
#' @keywords internal
quadrant_counts <- function(x,
                            y) {
    c(sum(x > 0 & y > 0), sum(x < 0 & y > 0), sum(x < 0 & y < 0), sum(x > 0 & y < 0))
}
