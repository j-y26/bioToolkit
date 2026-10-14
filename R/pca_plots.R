# Sample-level views of bulk data: PCA scores and variance explained by
# covariates (roadmap items B2, B4)
#
# Design
# ------
# * `plot_pca()` (B2) takes a `prcomp` result or a normalized features x
#   samples matrix (e.g. DESeq2 vst / rlog, edgeR logCPM). From a matrix it
#   keeps the `n.top` most variable features (500, as DESeq2::plotPCA()),
#   since thousands of flat, noisy features dilute the structure. A matrix of
#   raw counts is flagged: its PCA is driven by library size and the few most
#   highly expressed genes.
# * Samples are filled circles (or other filled shapes for `shape.by`) with a
#   thin dark outline, as in the composition scatter, so that light colours
#   and overlapping samples stay visible. Axis titles give each PC's share of
#   the variance (of the features used). Square panel (F2); the two axes keep
#   their own scales, as in most PCA plots: the % in the titles say how much
#   each axis counts.
# * Ellipses (`ellipse = TRUE`) are drawn only for groups of >= 4 samples: an
#   ellipse fitted to 2-3 points is mostly noise. Others are named in a
#   message.
# * `scree = TRUE` adds a small bar panel of the % variance of the first PCs
#   to the right, the PCs drawn in dark. (An inset scree, in the roadmap,
#   was left out: it covers samples in a corner of the panel.)
# * `plot_variance_explained()` (B4) has two views:
#     "genes": the share of each feature's variance explained by each
#              covariate, one violin (with a box) per covariate, ordered by
#              median, Residuals last in grey; or stacked bars for a few
#              `features`. From a variancePartition result (fractions from
#              linear mixed models, Hoffman & Schadt 2016, BMC Bioinformatics
#              17:483), or, from a matrix, a built-in approximation with
#              ordinary linear models (`variance_fractions()`).
#     "pcs":   R^2 of each top PC regressed on each covariate, as a grid of
#              tiles under bars of each PC's % variance: shows at once which
#              covariates drive the main axes of variation (e.g. a batch on
#              PC1).
# * The built-in "genes" fractions split each feature's R^2 (all covariates
#   in one linear model) between the covariates by the LMG / Shapley method:
#   a covariate's share is its increase in R^2 when added to the model,
#   averaged over all orders in which the covariates can enter (Lindeman,
#   Merenda & Gold 1980; Groemping 2006, R package relaimpo). Residuals =
#   1 - R^2, so each feature's fractions sum to 1, as variancePartition's,
#   whatever the order of the covariates, and variance shared by correlated
#   covariates is split between them. (First tried: the variance explained
#   beyond all other covariates, a semi-partial R^2; with an unbalanced
#   design its fractions summed to more than 1 for some genes, so they could
#   not be stacked.) Fixed effects only: with few samples a covariate with
#   many levels (e.g. individual) overfits; variancePartition's random
#   effects are better there.


#' @title PCA score plot of samples
#'
#' @description Plots the samples on two principal components, coloured,
#'              shaped and labelled by sample covariates, with the % of
#'              variance of each PC in the axis titles, optional ellipses per
#'              group and an optional scree panel. Takes a \code{prcomp}
#'              result, or a normalized features x samples matrix, of which
#'              the most variable features are used. Square panel of a fixed
#'              size; save it with \code{\link{save_plot}}; the scores and
#'              covariates behind it are returned by
#'              \code{\link{export_source_data}}.
#'
#' @param x A \code{prcomp} result (samples as rows of `x$x`), or a
#'          normalized features x samples matrix or data frame (e.g.
#'          \code{SummarizedExperiment::assay(DESeq2::vst(dds))} or
#'          \code{edgeR::cpm(y, log = TRUE)}), with sample names as column
#'          names. Raw counts should not be used: their PCA reflects library
#'          size and the few most highly expressed genes
#'
#' @param metadata Data frame of sample covariates, with sample names as row
#'                 names (in any order; extra rows are ignored). Needed for
#'                 `color.by`, `shape.by` and `label.by` other than TRUE.
#'                 Default is NULL
#'
#' @param color.by Column of `metadata` to colour samples by: discrete
#'                 (character, factor, logical) or continuous (numeric).
#'                 Default is NULL (one colour)
#'
#' @param shape.by Column of `metadata` (discrete, at most 5 levels) to give
#'                 samples different (filled) shapes. Default is NULL
#'
#' @param label.by Labels of the samples (needs ggrepel): TRUE for the sample
#'                 names, a column of `metadata`, or NULL (default) for none
#'
#' @param dims Two PCs to plot, default `c(1, 2)`
#'
#' @param n.top Matrix input: number of most variable features used, default
#'              500 (as \code{DESeq2::plotPCA()}); NULL or Inf for all
#'
#' @param center,scale Matrix input: whether features are centred (default
#'                     TRUE) and scaled to unit variance (default FALSE, so
#'                     that features keep their own variance, as is usual
#'                     for log-scale expression) before the PCA
#'
#' @param ellipse Whether to draw an ellipse around the samples of each
#'                `color.by` group (discrete), for groups of at least 4
#'                samples. Default is FALSE
#'
#' @param ellipse.level Confidence level of the ellipses, default 0.95
#'
#' @param ellipse.type "t" (default; multivariate t, robust to a stray
#'                     sample) or "norm" (multivariate normal), as in
#'                     \code{ggplot2::stat_ellipse()}
#'
#' @param scree Whether to add a panel of the % variance of the first PCs to
#'              the right, the plotted PCs in dark. Default is FALSE
#'
#' @param scree.n Number of PCs in the scree panel, default 10
#'
#' @param cols Colours of `color.by`: for a discrete covariate NULL
#'             (default; "scanpy-20"), a qualitative palette name, or colours
#'             named by level; for a continuous one NULL (default;
#'             "seq-blue"), a palette name or colours from low to high. With
#'             no `color.by`, one colour (default "grey40")
#'
#' @param pt.size Size of the points, default 3
#'
#' @param label.size Size of the labels in mm. If NULL (default), the small
#'                   text size of the theme (9 pt)
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
#' The % variance of a PC is its share of the total variance of the data the
#' PCA was run on, i.e. of the `n.top` features used for a matrix. Messages
#' report the number of features and samples used.
#'
#' The source data (\code{\link{export_source_data}}) hold one row per
#' sample: the sample name, the scores on the first PCs (up to 10) and the
#' covariates; with `scree = TRUE`, a second table holds the % variance per
#' PC.
#'
#' @return A ggplot object, or a patchwork with `scree = TRUE`
#'
#' @examples
#' set.seed(1)
#' m <- matrix(rnorm(2000 * 8, mean = 6), nrow = 2000,
#'             dimnames = list(paste0("Gene", 1:2000), paste0("S", 1:8)))
#' m[1:200, 5:8] <- m[1:200, 5:8] + 2    # a treatment effect
#' meta <- data.frame(treatment = rep(c("Control", "Treated"), each = 4),
#'                    batch = rep(c("A", "B"), 4),
#'                    row.names = colnames(m))
#' plot_pca(m, meta, color.by = "treatment", shape.by = "batch")
#' plot_pca(m, meta, color.by = "treatment", label.by = TRUE, scree = TRUE)
#' \dontrun{
#' vsd <- DESeq2::vst(dds)
#' p <- plot_pca(SummarizedExperiment::assay(vsd), as.data.frame(SummarizedExperiment::colData(vsd)),
#'               color.by = "condition", ellipse = TRUE)
#' save_plot(p, "pca.pdf")
#' }
#'
#' @seealso \code{\link{plot_variance_explained}}, \code{\link{save_plot}},
#'          \code{\link{export_source_data}}
#'
#' @import ggplot2
#' @export
#'
plot_pca <- function(
    x,
    metadata = NULL,
    color.by = NULL,
    shape.by = NULL,
    label.by = NULL,
    dims = c(1, 2),
    n.top = 500,
    center = TRUE,
    scale = FALSE,
    ellipse = FALSE,
    ellipse.level = 0.95,
    ellipse.type = c("t", "norm"),
    scree = FALSE,
    scree.n = 10,
    cols = NULL,
    pt.size = 3,
    label.size = NULL,
    panel.width = 3,
    panel.height = panel.width,
    legend.position = "right",
    bg = "transparent",
    seed = 42
) {
    ellipse.type <- match.arg(ellipse.type)
    if (!(is.numeric(dims) && length(dims) == 2 && !anyNA(dims) && all(dims >= 1) && dims[1] != dims[2])) {
        stop("`dims` must be two different PC numbers", call. = FALSE)
    }
    dims <- as.integer(dims)
    pc <- pca_input(x, n.top, center, scale)
    n_pc <- ncol(pc$scores)
    if (max(dims) > n_pc) stop("Only ", n_pc, " PCs are available", call. = FALSE)
    samples <- rownames(pc$scores)
    meta <- pca_metadata(metadata, samples)
    for (v in c("color.by", "shape.by")) {
        val <- get(v)
        if (!is.null(val) && !(is.character(val) && length(val) == 1 && val %in% names(meta))) {
            stop("`", v, "` must be a column of `metadata`", call. = FALSE)
        }
    }

    # One row per sample: its scores on the first PCs, then the covariates
    keep_pcs <- seq_len(min(n_pc, max(10, dims)))
    df <- data.frame(sample = samples, pc$scores[, keep_pcs, drop = FALSE], check.names = FALSE,
                     stringsAsFactors = FALSE)
    df <- cbind(df, meta[setdiff(names(meta), names(df))])
    rownames(df) <- NULL
    pc_x <- colnames(pc$scores)[dims[1]]
    pc_y <- colnames(pc$scores)[dims[2]]
    df$.x <- df[[pc_x]]
    df$.y <- df[[pc_y]]

    sz <- theme_sizes()
    lw <- pt_to_linewidth(sz$line)
    if (is.null(label.size)) label.size <- pt_to_size(sz$small)
    continuous <- !is.null(color.by) && is.numeric(df[[color.by]])

    # Fill: the colour covariate, or one colour
    if (is.null(color.by)) {
        df$.fill <- factor("all")
        fill_scale <- scale_fill_manual(values = c(all = if (is.null(cols)) "grey40" else cols[1]), guide = "none")
    } else if (continuous) {
        df$.fill <- df[[color.by]]
        fill_scale <- scale_fill_gradientn(
            colours = resolve_continuous_palette(cols, FALSE, default = "seq-blue"),
            name = color.by,
            guide = guide_colourbar(theme = theme(legend.key.height = unit(0.9, "in"),
                                                  legend.key.width = unit(0.12, "in"),
                                                  legend.frame = element_rect(colour = "black", linewidth = 0.3),
                                                  legend.ticks = element_line(colour = "black", linewidth = 0.3),
                                                  legend.ticks.length = unit(0.03, "in")))
        )
    } else {
        df$.fill <- as_group_factor(df[[color.by]])
        gc <- group_colours(cols, df$.fill)
        fill_scale <- scale_fill_manual(values = gc$cols, name = color.by, drop = TRUE)
    }
    # Shapes: filled ones (21-25), so that the fill colour and outline stay
    shape_scale <- NULL
    if (!is.null(shape.by)) {
        df$.shape <- as_group_factor(df[[shape.by]])
        n_shape <- nlevels(droplevels(df$.shape))
        if (n_shape > 5) {
            stop("`shape.by` has ", n_shape, " levels; at most 5 shapes can be told apart. ",
                 "Use it as `color.by` instead", call. = FALSE)
        }
        shape_scale <- scale_shape_manual(values = c(21, 24, 22, 23, 25), name = shape.by, drop = TRUE)
    }

    p <- ggplot(df, aes(x = .data$.x, y = .data$.y))
    if (ellipse) {
        if (is.null(color.by) || continuous) {
            message("`ellipse` needs a discrete `color.by`; no ellipses drawn")
        } else {
            n_per <- table(droplevels(df$.fill))
            small <- names(n_per)[n_per < 4]
            if (length(small)) {
                message("No ellipse for group(s) with fewer than 4 samples: ", paste(small, collapse = ", "))
            }
            de <- df[df$.fill %in% names(n_per)[n_per >= 4], , drop = FALSE]
            if (nrow(de)) {
                p <- p + stat_ellipse(
                    data = de, aes(fill = .data$.fill, colour = .data$.fill, group = .data$.fill),
                    geom = "polygon", type = ellipse.type, level = ellipse.level, alpha = 0.12,
                    linewidth = lw, show.legend = FALSE
                ) +
                    scale_colour_manual(values = gc$cols, guide = "none")
            }
        }
    }
    point_aes <- if (is.null(shape.by)) aes(fill = .data$.fill) else aes(fill = .data$.fill, shape = .data$.shape)
    # Filled circles with a thin dark outline; a fixed shape only when shapes
    # are not mapped (a fixed one would override the mapping)
    pt_args <- list(mapping = point_aes, size = pt.size, colour = "grey15", stroke = pt_to_stroke(0.5))
    if (is.null(shape.by)) pt_args$shape <- 21
    p <- p + do.call(geom_point, pt_args)
    p <- p + fill_scale + shape_scale
    if (!is.null(shape.by) && !continuous) {
        # Fill legend keys must be a filled shape to show their colour
        p <- p + guides(fill = guide_legend(override.aes = list(shape = 21), order = 1),
                        shape = guide_legend(override.aes = list(fill = "grey60"), order = 2))
    }

    if (!is.null(label.by)) {
        check_installed("ggrepel", reason = "to label samples")
        df$.label <- if (isTRUE(label.by)) {
            df$sample
        } else if (is.character(label.by) && length(label.by) == 1 && label.by %in% names(meta)) {
            as.character(df[[label.by]])
        } else {
            stop("`label.by` must be TRUE (sample names) or a column of `metadata`", call. = FALSE)
        }
        p <- p + ggrepel::geom_text_repel(
            data = df, aes(label = .data$.label), size = label.size, colour = "black",
            box.padding = 0.35, point.padding = 0.25, min.segment.length = 0.2,
            segment.colour = "grey50", segment.size = lw, max.overlaps = Inf, seed = seed, max.time = 2
        )
    }

    pct <- pc$variance * 100
    axis_title <- function(i) sprintf("%s (%.1f%%)", colnames(pc$scores)[i], pct[i])
    p <- p +
        scale_x_continuous(expand = expansion(mult = 0.08)) +
        scale_y_continuous(expand = expansion(mult = 0.08)) +
        labs(x = axis_title(dims[1]), y = axis_title(dims[2])) +
        theme_bioToolkit("axis", bg = bg) +
        panel_size_theme(panel.width, panel.height) +
        theme(legend.position = legend.position)
    src <- df[setdiff(names(df), c(".x", ".y", ".fill", ".shape", ".label"))]
    attr(p, "source_data") <- src
    if (!scree) return(p)

    # Scree: the % variance of the first PCs, the plotted ones in dark
    n_scree <- min(scree.n, n_pc)
    sd <- data.frame(pc = factor(seq_len(n_scree)), variance = pct[seq_len(n_scree)],
                     plotted = seq_len(n_scree) %in% dims)
    s <- ggplot(sd, aes(x = .data$pc, y = .data$variance, fill = .data$plotted)) +
        geom_col(width = 0.75) +
        scale_fill_manual(values = c(`TRUE` = "grey25", `FALSE` = "grey75"), guide = "none") +
        scale_x_discrete(expand = grid_expansion()) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                           labels = function(v) paste0(v, "%")) +
        labs(x = "PC", y = "Variance") +
        theme_bioToolkit("axis", bg = bg) +
        panel_size_theme(grid_panel_size(n_scree, 0.13), panel.height)
    out <- patchwork::wrap_plots(p, s, nrow = 1) +
        patchwork::plot_layout(guides = "collect") +
        patchwork::plot_annotation(theme = theme_bg(bg) + theme(legend.position = legend.position))
    out
}


#' @title Variance explained by sample covariates
#'
#' @description Shows how much of the variation in a dataset each sample
#'              covariate (condition, batch, sex, library size, ...)
#'              explains, in one of two views:
#'              \describe{
#'                \item{"genes"}{The share of each feature's variance explained
#'                      by each covariate: one violin (with a box of the
#'                      quartiles and median) per covariate, ordered by
#'                      median, Residuals last; or, for a few `features`,
#'                      stacked bars. From the result of
#'                      \code{variancePartition::fitExtractVarPartModel()}, or
#'                      computed here with linear models (see Details).}
#'                \item{"pcs"}{The R-squared of each of the top PCs regressed
#'                      on each covariate, as tiles under bars of each PC's \%
#'                      variance: which covariates drive the main axes of
#'                      variation.}
#'              }
#'              Fixed sizes; save it with \code{\link{save_plot}}; the values
#'              behind it are returned by \code{\link{export_source_data}}.
#'
#' @param x A normalized features x samples matrix or data frame (e.g. vst,
#'          logCPM), sample names as column names; or, for view "genes", a
#'          variancePartition result (features x components, with a
#'          `Residuals` column); or, for view "pcs", a \code{prcomp} result
#'
#' @param metadata Data frame of sample covariates, sample names as row
#'                 names. Not needed for a variancePartition result
#'
#' @param covariates Columns of `metadata` to use. If NULL (default), all
#'
#' @param view "genes" (default) or "pcs"
#'
#' @param n.top Number of most variable features used from a matrix. If NULL
#'              (default): all features for view "genes", 500 for view "pcs"
#'              (as \code{\link{plot_pca}})
#'
#' @param n.pcs View "pcs": number of PCs, default 10
#'
#' @param features View "genes": features to draw as stacked bars instead of
#'                 violins (e.g. a few genes of interest). Default is NULL
#'
#' @param cols Colours of the covariates: NULL (default; "scanpy-20"), a
#'             qualitative palette name, or colours named by covariate
#'             (Residuals are always grey). For view "pcs", the colours of
#'             the R-squared scale: NULL (default; "seq-blue"), a palette
#'             name or colours from low to high
#'
#' @param italic.genes Whether feature names (stacked bars) are gene names,
#'                     set in italics. Default is TRUE
#'
#' @param cell.width Width of one violin (view "genes") or one PC column
#'                   (view "pcs") in inches, default 0.4
#'
#' @param cell.height Height of one bar (stacked bars) or one covariate row
#'                    (view "pcs") in inches, default 0.25
#'
#' @param panel.width Length of the stacked bars in inches, default 3
#'
#' @param panel.height Height of the violins in inches, default 2
#'
#' @param legend.position Position of the legend (stacked bars and view
#'                        "pcs"), default "right"
#'
#' @param bg Background colour, default is "transparent"
#'
#' @details
#' \strong{Built-in fractions (view "genes" from a matrix).} For each feature,
#' a linear model of its values on all covariates is fitted (numeric
#' covariates as they are, others as factors), and its R-squared is split
#' between the covariates by the LMG method (Lindeman, Merenda & Gold 1980;
#' as \code{relaimpo::calc.relimp(type = "lmg")}, Groemping 2006): a
#' covariate's fraction is its increase in R-squared when added to the
#' model, averaged over all orders in which the covariates can be added.
#' Residuals is 1 - R-squared, so the fractions of a feature sum to 1, as
#' variancePartition's, and do not depend on the order of `covariates`.
#' Variance shared by correlated covariates (e.g. a batch partly confounded
#' with condition) is split between them; fully confounded covariates
#' cannot be separated and give an error. This is an approximation of
#' variancePartition: fixed effects only. With few samples, a covariate
#' with many levels is overfitted; use variancePartition, with such
#' covariates as random effects, for those designs. At most 8 covariates.
#'
#' \strong{View "pcs".} Each PC's scores are regressed on each covariate
#' separately (R-squared of a one-covariate linear model; for a factor, the
#' share of the PC's variance between its levels). The source data also hold
#' the F-test p-value and its Benjamini-Hochberg adjustment over all tiles;
#' with few samples they are rough.
#'
#' @return A ggplot object, or a patchwork for view "pcs"
#'
#' @examples
#' set.seed(1)
#' m <- matrix(rnorm(1000 * 12, mean = 6), nrow = 1000,
#'             dimnames = list(paste0("Gene", 1:1000), paste0("S", 1:12)))
#' meta <- data.frame(condition = rep(c("WT", "KO"), each = 6),
#'                    batch = rep(c("A", "B", "C"), 4),
#'                    row.names = colnames(m))
#' m[1:300, meta$condition == "KO"] <- m[1:300, meta$condition == "KO"] + 1.5
#' m[, meta$batch == "B"] <- m[, meta$batch == "B"] + rnorm(1000, sd = 0.5)
#' plot_variance_explained(m, meta)
#' plot_variance_explained(m, meta, features = c("Gene1", "Gene2", "Gene900"))
#' plot_variance_explained(m, meta, view = "pcs")
#'
#' @seealso \code{\link{plot_pca}}, \code{\link{save_plot}}
#'
#' @import ggplot2
#' @export
#'
plot_variance_explained <- function(
    x,
    metadata = NULL,
    covariates = NULL,
    view = c("genes", "pcs"),
    n.top = NULL,
    n.pcs = 10,
    features = NULL,
    cols = NULL,
    italic.genes = TRUE,
    cell.width = 0.4,
    cell.height = 0.25,
    panel.width = 3,
    panel.height = 2,
    legend.position = "right",
    bg = "transparent"
) {
    view <- match.arg(view)
    sz <- theme_sizes()
    if (view == "pcs") {
        return(variance_pcs_plot(x, metadata, covariates, n.top, n.pcs, cols, cell.width, cell.height,
                                 legend.position, bg, sz))
    }

    # Fractions: features x components, Residuals last
    is_vp <- inherits(x, "varPartResults") || (is.data.frame(x) && "Residuals" %in% names(x))
    if (is_vp) {
        message("Reading variance fractions (a variancePartition result)")
        fr <- as.data.frame(x)
        fr <- fr[c(setdiff(names(fr), "Residuals"), "Residuals")]
        if (!is.null(covariates)) {
            missing <- setdiff(covariates, names(fr))
            if (length(missing)) stop("Not in the result: ", paste(missing, collapse = ", "), call. = FALSE)
            fr <- fr[c(covariates, "Residuals")]
        }
    } else {
        m <- as_feature_matrix(x)
        meta <- pca_metadata(metadata, colnames(m), required = TRUE)
        if (is.null(covariates)) covariates <- names(meta)
        m <- top_variable(m, n.top)
        fr <- variance_fractions(m, meta, covariates)
    }
    comps <- setdiff(names(fr), "Residuals")
    if (length(comps) == 0) stop("No covariates to show", call. = FALSE)
    # Covariates by decreasing median fraction; Residuals last, in grey
    med <- vapply(fr[comps], stats::median, numeric(1), na.rm = TRUE)
    lv <- c(comps[order(-med)], "Residuals")
    comp_cols <- group_colours(cols, factor(comps, levels = comps))$cols
    comp_cols <- c(comp_cols[comps], Residuals = "grey80")

    if (!is.null(features)) {
        return(variance_bars_plot(fr, features, lv, comp_cols, italic.genes, cell.height, panel.width,
                                  legend.position, bg, sz))
    }

    long <- data.frame(
        feature = rep(rownames(fr), times = length(lv)),
        group = factor(rep(lv, each = nrow(fr)), levels = lv),
        value = unlist(fr[lv], use.names = FALSE),
        stringsAsFactors = FALSE
    )
    long <- long[!is.na(long$value), , drop = FALSE]
    long$drawn <- TRUE
    # The violins and boxes of VlnPlot_bioToolkit(), without points
    p <- violin_panel(
        d = long,
        pts = long[0, , drop = FALSE],
        group_cols = comp_cols,
        flip = FALSE,
        box = TRUE,
        box.width = 0.2,
        pt.size = 0,
        alpha = 0,
        adjust = 1,
        layer_fun = identity,
        seed = 1,
        annotate = FALSE
    ) +
        scale_y_continuous(labels = scales::label_percent(drop0trailing = TRUE), breaks = seq(0, 1, 0.25),
                           expand = expansion(mult = c(0.02, 0.02))) +
        coord_cartesian(ylim = c(0, 1), clip = "off") +
        labs(x = NULL, y = "Variance explained") +
        group_axis_theme(FALSE, sz) +
        theme(legend.position = "none") +
        distribution_size_theme(long, FALSE, cell.width, NULL, NULL, panel.height) +
        theme_bg(bg)
    attr(p, "source_data") <- variance_source(fr)
    p
}


#' Stacked bars of the variance fractions of a few features
#'
#' @keywords internal
variance_bars_plot <- function(fr,
                               features,
                               lv,
                               comp_cols,
                               italic.genes,
                               cell.height,
                               panel.width,
                               legend.position,
                               bg,
                               sz) {
    features <- as.character(features)
    missing <- setdiff(features, rownames(fr))
    if (length(missing)) {
        near <- rownames(fr)[tolower(rownames(fr)) %in% tolower(missing)]
        warning("`features` not found: ", paste(missing, collapse = ", "),
                if (length(near)) paste0(" (matching only in case: ", paste(unique(near), collapse = ", "), ")"),
                call. = FALSE)
        features <- intersect(features, rownames(fr))
    }
    if (length(features) == 0) stop("None of `features` is in the data", call. = FALSE)
    sub <- fr[features, lv, drop = FALSE]
    long <- data.frame(
        feature = factor(rep(features, times = length(lv)), levels = rev(features)),
        component = factor(rep(lv, each = length(features)), levels = rev(lv)),
        value = unlist(sub, use.names = FALSE)
    )
    p <- ggplot(long, aes(x = .data$value, y = .data$feature, fill = .data$component)) +
        geom_col(width = 0.8, colour = NA) +
        scale_fill_manual(values = comp_cols, name = NULL, breaks = lv) +
        scale_x_continuous(labels = scales::label_percent(drop0trailing = TRUE), expand = expansion(0),
                           limits = c(0, max(1, rowSums(sub, na.rm = TRUE)))) +
        scale_y_discrete(expand = grid_expansion(), labels = function(v) italic_labels(v, italic.genes)) +
        labs(x = "Variance explained", y = NULL) +
        theme_bioToolkit("axis", bg = bg) +
        panel_size_theme(panel.width, grid_panel_size(length(features), cell.height)) +
        theme(legend.position = legend.position)
    attr(p, "source_data") <- variance_source(sub)
    p
}


#' The "pcs" view of `plot_variance_explained()`
#'
#' @keywords internal
variance_pcs_plot <- function(x,
                              metadata,
                              covariates,
                              n.top,
                              n.pcs,
                              cols,
                              cell.width,
                              cell.height,
                              legend.position,
                              bg,
                              sz) {
    if (is.null(n.top)) n.top <- 500
    pc <- pca_input(x, n.top, center = TRUE, scale = FALSE)
    meta <- pca_metadata(metadata, rownames(pc$scores), required = TRUE)
    if (is.null(covariates)) covariates <- names(meta)
    check_covariates(meta, covariates)
    n_pc <- min(n.pcs, ncol(pc$scores))
    pcs <- colnames(pc$scores)[seq_len(n_pc)]

    # R^2 and F-test p of PC ~ covariate, one model per tile
    tiles <- do.call(rbind, lapply(covariates, function(cv) {
        v <- covariate_column(meta[[cv]])
        do.call(rbind, lapply(pcs, function(k) {
            ok <- !is.na(v)
            fit <- stats::lm(pc$scores[ok, k] ~ v[ok])
            s <- summary(fit)
            f <- s$fstatistic
            p <- if (is.null(f)) NA_real_ else stats::pf(f[1], f[2], f[3], lower.tail = FALSE)
            data.frame(covariate = cv, pc = k, r2 = s$r.squared, p = unname(p))
        }))
    }))
    tiles$padj <- stats::p.adjust(tiles$p, "BH")
    tiles$covariate <- factor(tiles$covariate, levels = rev(covariates))
    tiles$pc <- factor(tiles$pc, levels = pcs)
    tiles$variance <- pc$variance[match(as.character(tiles$pc), colnames(pc$scores))] * 100
    # Default: "seq-blue" from near white, so that tiles of R^2 near 0 recede
    fill_cols <- if (is.null(cols)) c("grey97", bt_palette("seq-blue")) else resolve_continuous_palette(cols, FALSE)
    # Values in white on dark tiles, dark on mid ones, light grey below 0.1
    tiles$label <- formatC(tiles$r2, format = "f", digits = 2)
    tiles$text <- cut(tiles$r2, c(-Inf, 0.1, 0.55, Inf), labels = c("low", "mid", "high"))

    width <- grid_panel_size(n_pc, cell.width)
    text_size <- pt_to_size(sz$small * 0.85)
    heat <- ggplot(tiles, aes(x = .data$pc, y = .data$covariate, fill = .data$r2)) +
        geom_tile(colour = if (is.na(bg) || identical(bg, "transparent")) "white" else bg,
                  linewidth = pt_to_linewidth(1)) +
        geom_text(aes(label = .data$label, colour = .data$text), size = text_size, show.legend = FALSE) +
        scale_colour_manual(values = c(low = "grey60", mid = "grey15", high = "white"), guide = "none") +
        scale_fill_gradientn(
            colours = fill_cols, limits = c(0, 1), name = quote(R^2),
            guide = guide_colourbar(theme = theme(legend.key.height = unit(0.9, "in"),
                                                  legend.key.width = unit(0.12, "in"),
                                                  legend.frame = element_rect(colour = "black", linewidth = 0.3),
                                                  legend.ticks = element_line(colour = "black", linewidth = 0.3),
                                                  legend.ticks.length = unit(0.03, "in")))
        ) +
        scale_x_discrete(expand = expansion(0)) +
        scale_y_discrete(expand = expansion(0)) +
        theme_bioToolkit("grid", bg = bg) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(),
              axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
              legend.position = legend.position) +
        panel_size_theme(width, grid_panel_size(length(covariates), cell.height))
    attr(heat, "source_data") <- tiles[c("pc", "variance", "covariate", "r2", "p", "padj")]

    # Bars of each PC's % variance above, sharing the columns
    vd <- data.frame(pc = factor(pcs, levels = pcs), variance = pc$variance[seq_len(n_pc)] * 100)
    bars <- ggplot(vd, aes(x = .data$pc, y = .data$variance)) +
        geom_col(width = 0.7, fill = "grey55") +
        geom_text(aes(label = sprintf("%.0f%%", .data$variance)), vjust = -0.4, size = text_size,
                  colour = "grey20") +
        scale_x_discrete(expand = expansion(0)) +
        # The name of the bars as the one tick label of the y-axis, so that it
        # lines up with the covariate names below (an axis title would sit
        # left of them)
        scale_y_continuous(expand = expansion(mult = c(0, 0.3)), breaks = max(vd$variance) / 2,
                           labels = "% variance") +
        labs(x = NULL, y = NULL) +
        theme_bioToolkit("axis", bg = bg) +
        theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line.y = element_blank(),
              axis.text.y = element_text(size = sz$small, colour = "black")) +
        panel_size_theme(width, 0.5)
    patchwork::wrap_plots(bars, heat, ncol = 1) +
        patchwork::plot_layout(guides = "collect") +
        patchwork::plot_annotation(theme = theme_bg(bg) + theme(legend.position = legend.position))
}


#' Variance fractions of each feature, from linear models
#'
#' For each feature, the share of its variance explained by each covariate
#' (the LMG / Shapley decomposition of R^2: the increase in R^2 when the
#' covariate is added, averaged over all orders in which the covariates can
#' enter the model) and the residual fraction (1 - R^2). All features are
#' fitted at once, one QR decomposition per subset of covariates.
#'
#' @param m Features x samples matrix
#'
#' @param meta Data frame of covariates, rows in the order of the columns of
#'             `m`
#'
#' @param covariates Columns of `meta`
#'
#' @return A data frame, features x (covariates, Residuals); each row sums
#'         to 1
#'
#' @references Lindeman, Merenda & Gold (1980) Introduction to Bivariate and
#'             Multivariate Analysis. Groemping (2006) Relative importance
#'             for linear regression in R: the package relaimpo. J Stat Softw
#'             17(1)
#'
#' @keywords internal
variance_fractions <- function(m,
                               meta,
                               covariates) {
    check_covariates(meta, covariates)
    k <- length(covariates)
    if (k > 8) stop("At most 8 covariates (2^k models are fitted)", call. = FALSE)
    d <- data.frame(lapply(meta[covariates], covariate_column), check.names = FALSE)
    ok <- stats::complete.cases(d)
    if (!all(ok)) {
        message(sum(!ok), " sample(s) with missing covariates are left out")
        d <- d[ok, , drop = FALSE]
        m <- m[, ok, drop = FALSE]
    }
    # Features that do not vary have no variance to explain
    tss <- rowSums((m - rowMeans(m))^2)
    flat <- !(tss > 0)
    if (any(flat)) {
        message(sum(flat), " feature(s) without variance are left out")
        m <- m[!flat, , drop = FALSE]
        tss <- tss[!flat]
    }
    model_matrix <- function(cv) {
        if (length(cv) == 0) return(matrix(1, nrow(d), 1))
        stats::model.matrix(stats::reformulate(paste0("`", cv, "`")), d)
    }
    X <- model_matrix(covariates)
    rank_full <- qr(X)$rank
    if (rank_full < ncol(X)) {
        stop("The covariates are confounded (", paste(covariates, collapse = ", "),
             "): some cannot be separated from the others. Drop or merge them", call. = FALSE)
    }
    if (nrow(d) - rank_full < 1) {
        stop("No residual degrees of freedom: ", nrow(d), " samples for ", rank_full,
             " model coefficients. Use fewer covariates", call. = FALSE)
    }
    # R^2 of every subset of covariates: subset i holds the covariates whose
    # bits are set in i - 1 (subset 1 is the intercept only, R^2 = 0)
    bits <- 2^(seq_len(k) - 1)
    subsets <- lapply(seq_len(2^k) - 1, function(b) covariates[bitwAnd(b, bits) > 0])
    r2 <- vapply(subsets, function(cv) {
        1 - colSums(qr.resid(qr(model_matrix(cv)), t(m))^2) / tss
    }, numeric(nrow(m)))
    r2 <- matrix(r2, nrow = nrow(m))
    # Share of covariate j: the increase in R^2 from adding j to a subset S
    # of the others, weighted by |S|! (k - |S| - 1)! / k!, the fraction of
    # the k! orders in which exactly the covariates of S come before j
    sizes <- vapply(subsets, length, integer(1))
    out <- vapply(seq_len(k), function(j) {
        without <- which(bitwAnd(seq_len(2^k) - 1, bits[j]) == 0)
        w <- factorial(sizes[without]) * factorial(k - sizes[without] - 1) / factorial(k)
        drop((r2[, without + bits[j], drop = FALSE] - r2[, without, drop = FALSE]) %*% w)
    }, numeric(nrow(m)))
    out <- matrix(out, nrow = nrow(m), dimnames = list(rownames(m), covariates))
    data.frame(out, Residuals = 1 - r2[, 2^k], check.names = FALSE)
}


#' Check covariates for the variance models
#'
#' @keywords internal
check_covariates <- function(meta,
                             covariates) {
    missing <- setdiff(covariates, names(meta))
    if (length(missing)) stop("Not columns of `metadata`: ", paste(missing, collapse = ", "), call. = FALSE)
    for (cv in covariates) {
        v <- meta[[cv]]
        n <- length(unique(v[!is.na(v)]))
        if (n < 2) stop("Covariate '", cv, "' has a single value", call. = FALSE)
        if (!is.numeric(v) && n == sum(!is.na(v))) {
            stop("Covariate '", cv, "' has a different level for every sample; it explains everything",
                 call. = FALSE)
        }
    }
    invisible(TRUE)
}


#' A covariate as a model term: numbers as they are, the rest as factors
#'
#' @keywords internal
covariate_column <- function(v) {
    if (is.numeric(v)) v else droplevels(as.factor(v))
}


#' Source data of the "genes" view: one row per feature
#'
#' @keywords internal
variance_source <- function(fr) {
    data.frame(feature = rownames(fr), fr, check.names = FALSE, row.names = NULL)
}


#' PCA of a matrix or a `prcomp` result
#'
#' @param x A `prcomp` result, or a features x samples matrix / data frame
#'
#' @param n.top Most variable features to keep (NULL or Inf for all)
#'
#' @param center,scale Passed to `prcomp()`
#'
#' @return A list with `scores` (samples x PCs, columns "PC1", ...) and
#'         `variance` (share of the total variance per PC)
#'
#' @keywords internal
pca_input <- function(x,
                      n.top,
                      center,
                      scale) {
    if (inherits(x, "prcomp")) {
        scores <- x$x
        if (is.null(rownames(scores))) stop("The PCA scores have no sample names", call. = FALSE)
        colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
        return(list(scores = scores, variance = x$sdev^2 / sum(x$sdev^2)))
    }
    m <- as_feature_matrix(x)
    m <- top_variable(m, n.top)
    # Features without variance cannot be scaled, and add nothing
    v <- apply(m, 1, stats::var)
    m <- m[v > 0, , drop = FALSE]
    if (nrow(m) < 2) stop("Fewer than 2 features vary across samples", call. = FALSE)
    message("PCA of ", nrow(m), " features x ", ncol(m), " samples")
    pr <- stats::prcomp(t(m), center = center, scale. = scale)
    scores <- pr$x
    colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
    list(scores = scores, variance = pr$sdev^2 / sum(pr$sdev^2))
}


#' A features x samples matrix from a matrix or data frame
#'
#' @keywords internal
as_feature_matrix <- function(x) {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (!(is.matrix(x) && is.numeric(x))) {
        stop("`x` must be a numeric features x samples matrix (or a prcomp result)", call. = FALSE)
    }
    if (is.null(colnames(x))) stop("`x` has no sample (column) names", call. = FALSE)
    if (anyDuplicated(colnames(x))) stop("Sample (column) names of `x` are not unique", call. = FALSE)
    if (ncol(x) < 3) stop("At least 3 samples are needed", call. = FALSE)
    if (anyNA(x)) {
        keep <- stats::complete.cases(x)
        message(sum(!keep), " feature(s) with missing values are left out")
        x <- x[keep, , drop = FALSE]
    }
    # Raw counts: large non-negative integers. Their variance is dominated by
    # library size and the most highly expressed genes
    if (all(x >= 0) && max(x) > 1000 && all(x == round(x))) {
        warning("`x` looks like raw counts; use normalized values on a log scale (e.g. DESeq2::vst(), ",
                "edgeR::cpm(log = TRUE))", call. = FALSE)
    }
    x
}


#' The most variable features of a matrix
#'
#' @keywords internal
top_variable <- function(m,
                         n.top) {
    if (is.null(n.top) || is.infinite(n.top) || n.top >= nrow(m)) return(m)
    if (!(is.numeric(n.top) && length(n.top) == 1 && n.top >= 2)) {
        stop("`n.top` must be a single number >= 2, or NULL for all", call. = FALSE)
    }
    v <- apply(m, 1, stats::var)
    m[order(-v)[seq_len(n.top)], , drop = FALSE]
}


#' Sample covariates in the order of the samples
#'
#' @param metadata Data frame with sample names as row names, or NULL
#'
#' @param samples Sample names
#'
#' @param required Whether `metadata` must be given
#'
#' @return A data frame with one row per sample (no columns when `metadata`
#'         is NULL)
#'
#' @keywords internal
pca_metadata <- function(metadata,
                         samples,
                         required = FALSE) {
    if (is.null(metadata)) {
        if (required) stop("`metadata` (sample covariates) is needed", call. = FALSE)
        return(data.frame(row.names = samples))
    }
    metadata <- as.data.frame(metadata)
    missing <- setdiff(samples, rownames(metadata))
    if (length(missing)) {
        stop("Sample(s) not in the row names of `metadata`: ",
             paste(utils::head(missing, 5), collapse = ", "), if (length(missing) > 5) ", ...", call. = FALSE)
    }
    metadata[samples, , drop = FALSE]
}


#' A discrete covariate as a factor, keeping the order of factor levels
#'
#' @keywords internal
as_group_factor <- function(v) {
    if (is.factor(v)) droplevels(v) else factor(v, levels = unique(sort(v)))
}
