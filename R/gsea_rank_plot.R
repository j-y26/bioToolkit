# GSEA running enrichment score plot (roadmap item G3, addition 2026-10-13)
#
# Design
# ------
# `plot_gsea_rank()` replaces enrichplot::gseaplot2() as used in the user's
# own analyses. gseaplot2() returns an aplot/gglist object whose panels are
# hard to restyle, and needs enrichplot (not a dependency here). The running
# enrichment score is simple to compute, so it is computed here:
#   * genes are ranked by the metric, largest first;
#   * walking down the list, the score rises at each gene of the set by
#     |metric|^p / (sum of |metric|^p over the set's genes) and falls at each
#     other gene by 1 / (number of other genes);
#   * the enrichment score (ES) is the largest deviation from 0.
# This is the weighted Kolmogorov-Smirnov statistic of Subramanian et al.
# (2005) PNAS 102:15545, with p = 1 as in GSEA, fgsea (`gseaParam = 1`) and
# clusterProfiler (`exponent = 1`). The ES computed is checked against the
# one reported, to catch a ranked list that is not the one tested.
#
# Three stacked panels of one width, touching (no margins between them, one
# continuous left axis line), combined with patchwork, compact by default
# (1.3 / 0.32 / 0.55 in):
#   1. the running score, line coloured by the sign of the NES;
#   2. the positions of the set's genes (ticks) above a band coloured by the
#      ranking metric (blue - white - magenta). enrichplot::gseaplot2() draws
#      a band there too, but its colours only number 8 stretches of the list
#      holding equal shares of the set's genes, so they encode no value; here
#      the band shows where the list turns from up- to down-ranked genes, as
#      in the GSEA desktop plot;
#   3. the ranking metric (log2FC, signed -log10 p, t...), the same for every
#      gene set of a comparison, so drawn quietly in grey; `show.metric =
#      FALSE` leaves it out.
# NES and FDR (the adjusted p-value) sit in the corner the curve leaves
# empty: top right for a positive score, bottom left for a negative one. The title is the gene
# set's name; the subtitle its ID (when it differs) and the comparison, in
# italics; plain ggplot2 title and subtitle, so no ggtext.


#' @title GSEA running enrichment score plot
#'
#' @description Draws the running enrichment score of one gene set along the
#'              ranked gene list, the positions of the set's genes above a
#'              band coloured by the ranking metric (where genes rank up or
#'              down), and the ranking metric, with the NES and adjusted
#'              p-value (FDR) of the set. Takes a clusterProfiler \code{GSEA()} result directly, or
#'              fgsea (or any \code{\link{standardize_enrichment}}) results
#'              with the ranked list and the gene sets. Panels have a fixed
#'              size; save the plot with \code{\link{save_plot}}.
#'
#' @param gsea A clusterProfiler `gseaResult` (holds the ranked list and gene
#'             sets), or a GSEA results table (e.g. from
#'             \code{fgsea::fgsea()}) or a named list of them, as read by
#'             \code{\link{standardize_enrichment}}
#'
#' @param term The gene set: its ID (e.g. "GOBP_RESPONSE_TO_KETONE",
#'             "GO:1901654") or its name, original or cleaned (e.g. "Response
#'             to ketone"). When a name matches several gene sets, the most
#'             significant is drawn, with a message
#'
#' @param stats The ranked list: a numeric vector of the ranking metric
#'              (e.g. log2 fold change, signed -log10 p-value) named by gene.
#'              Needed unless `gsea` is a `gseaResult`
#'
#' @param gene.sets Gene sets: a named list of character vectors (names = the
#'                  IDs in `gsea`, e.g. `pathway` of fgsea). Needed unless
#'                  `gsea` is a `gseaResult`
#'
#' @param comparison The comparison: when `gsea` holds several, the one to
#'                   draw; otherwise only its name, written in the subtitle.
#'                   Default is NULL
#'
#' @param title Title. If NULL (default), the gene set's (cleaned) name
#'
#' @param show.id Whether the subtitle gives the gene set's ID when it differs
#'                from its name, default is TRUE
#'
#' @param cols Colours of the line and metric for a positive and a negative
#'             score: two colours (positive, negative), or named `up` /
#'             `down`. If NULL (default), the ends of "blue-magenta"
#'
#' @param metric.title Title of the ranked metric's axis, default is
#'                     "Ranked\\nmetric". Name the metric when it helps,
#'                     e.g. "log2 fold\\nchange"
#'
#' @param show.metric Whether to draw the ranked metric below the gene
#'                    positions, default is TRUE. It is the same for every
#'                    gene set of a comparison (the ranked list), so it can be
#'                    left out of all but one figure
#'
#' @param panel.width Width of the panels in inches, default is 3
#'
#' @param panel.heights Heights of the running score, gene position and
#'                      metric panels in inches, default is
#'                      `c(1.3, 0.32, 0.55)`
#'
#' @param bg Background colour, default is "transparent"
#'
#' @param ... Arguments to \code{\link{standardize_enrichment}} for a table
#'            (e.g. `comparison.col`)
#'
#' @return A patchwork object; its source data
#'         (\code{\link{export_source_data}}) hold one row per gene: `rank`,
#'         `gene`, `metric`, `running_es` and `in_set`
#'
#' @examples
#' set.seed(1)
#' genes <- paste0("Gene", 1:2000)
#' stats <- sort(stats::setNames(rnorm(2000), genes), decreasing = TRUE)
#' sets <- list(SET_UP = sample(genes[1:300], 40), SET_RANDOM = sample(genes, 40))
#' res <- data.frame(pathway = names(sets), pval = c(1e-4, 0.6), padj = c(2e-4, 0.6),
#'                   ES = NA, NES = c(2.1, 0.4), size = 40)
#' plot_gsea_rank(res, "SET_UP", stats = stats, gene.sets = sets, comparison = "KO vs WT")
#'
#' @seealso \code{\link{plot_enrichment_dots}},
#'          \code{\link{plot_enrichment_bars}},
#'          \code{\link{standardize_enrichment}}
#'
#' @import ggplot2
#' @export
#'
plot_gsea_rank <- function(
    gsea,
    term,
    stats = NULL,
    gene.sets = NULL,
    comparison = NULL,
    title = NULL,
    show.id = TRUE,
    cols = NULL,
    metric.title = "Ranked\nmetric",
    show.metric = TRUE,
    panel.width = 3,
    panel.heights = c(1.3, 0.32, 0.55),
    bg = "transparent",
    ...
) {
    if (!(is.character(term) && length(term) == 1 && !is.na(term))) {
        stop("`term` must be a single gene set ID or name", call. = FALSE)
    }
    if (!(is.numeric(panel.heights) && length(panel.heights) == 3 && all(panel.heights > 0))) {
        stop("`panel.heights` must be three positive numbers (inches)", call. = FALSE)
    }
    # The ranked list and gene sets come with a clusterProfiler result
    if (methods::is(gsea, "gseaResult")) {
        if (is.null(stats)) stats <- gsea@geneList
        if (is.null(gene.sets)) gene.sets <- gsea@geneSets
    }
    if (is.null(stats) || is.null(gene.sets)) {
        stop("Give the ranked list (`stats`) and the gene sets (`gene.sets`) the results were computed from",
             call. = FALSE)
    }
    if (!(is.numeric(stats) && !is.null(names(stats)))) stop("`stats` must be a numeric vector named by gene", call. = FALSE)
    if (!(is.list(gene.sets) && !is.null(names(gene.sets)))) {
        stop("`gene.sets` must be a named list of character vectors", call. = FALSE)
    }

    # The gene set's row: matched on its ID or name. For clusterProfiler the
    # standardized `term_id` is the description, so the ID (the name in
    # `gene.sets`) is read from the raw table
    res <- as_enrichment(gsea, ...)
    parts <- if (is.list(gsea) && !is.data.frame(gsea) && !isS4(gsea)) gsea else list(gsea)
    ids <- do.call(rbind, lapply(parts, function(g) {
        d <- enrichment_as_data_frame(g)
        id_col <- intersect(c("ID", "pathway", "NAME", "term_id"), names(d))[1]
        es_col <- intersect(c("enrichmentScore", "ES"), names(d))[1]
        data.frame(set_id = if (is.na(id_col)) NA_character_ else as.character(d[[id_col]]),
                   es = if (is.na(es_col)) NA_real_ else suppressWarnings(as.numeric(d[[es_col]])))
    }))
    # Rows of the standardized table follow those of the tables read
    res$set_id <- if (nrow(ids) == nrow(res) && !anyNA(ids$set_id)) ids$set_id else res$term_id
    res$es <- if (nrow(ids) == nrow(res)) ids$es else NA_real_
    if (!is.null(res$comparison) && nlevels(res$comparison) > 1) {
        if (is.null(comparison)) {
            stop("The results hold several comparisons (", paste(levels(res$comparison), collapse = ", "),
                 "); choose one with `comparison`", call. = FALSE)
        }
        if (!(comparison %in% levels(res$comparison))) stop("Comparison '", comparison, "' not found", call. = FALSE)
        res <- res[res$comparison == comparison, , drop = FALSE]
    } else if (is.null(comparison) && !is.null(res$comparison)) {
        comparison <- levels(droplevels(res$comparison))[1]
    }
    hits <- which(res$set_id == term | res$term_id == term | res$term == term)
    if (length(hits) == 0) stop("Gene set '", term, "' not found in the results", call. = FALSE)
    if (length(hits) > 1) {
        hits <- hits[order(res$padj[hits])]
        message(length(hits), " gene sets match '", term, "'; drawing the most significant, ", res$set_id[hits[1]])
    }
    info <- res[hits[1], , drop = FALSE]
    if (!(info$set_id %in% names(gene.sets))) {
        stop("Gene set '", info$set_id, "' is not in `gene.sets`", call. = FALSE)
    }

    # Running enrichment score
    run <- gsea_running_score(stats, gene.sets[[info$set_id]])
    peak <- which.max(abs(run$running_es))
    es_calc <- run$running_es[peak]
    if (!is.na(info$es) && abs(es_calc - info$es) > 0.01 * max(1, abs(info$es))) {
        warning("The enrichment score computed from `stats` (", signif(es_calc, 3), ") differs from the one ",
                "reported (", signif(info$es, 3), "): is `stats` the ranked list that was tested?", call. = FALSE)
    }
    nes <- if (is.na(info$nes)) es_calc else info$nes
    positive <- nes > 0

    # Colours: positive and negative scores
    ends <- bt_palette("blue-magenta")
    dir_cols <- c(up = ends[length(ends)], down = ends[1])
    if (!is.null(cols)) {
        if (is.null(names(cols))) names(cols) <- c("up", "down")[seq_along(cols)]
        dir_cols[names(cols)] <- cols
    }
    line_col <- if (positive) dir_cols[["up"]] else dir_cols[["down"]]

    sz <- theme_sizes()
    lw <- pt_to_linewidth(sz$line)
    n <- nrow(run)
    x_scale <- scale_x_continuous(limits = c(0.5, n + 0.5), expand = expansion(0))
    # Panels touch: no margin between them, and the left axis line runs
    # through all three
    stacked <- theme_bioToolkit("axis", bg = bg) +
        theme(
            axis.title.y = element_text(size = sz$small, lineheight = 0.9),
            axis.text = element_text(size = sz$small * 0.9)
        )
    no_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
                  axis.title.x = element_blank())

    # 1. Running score, NES and FDR in the corner the curve leaves empty
    stat_label <- paste0("NES = ", formatC(nes, format = "f", digits = 2), "\n",
                         "FDR = ", format(signif(info$padj, 2), scientific = info$padj < 1e-3))
    p1 <- ggplot(run, aes(x = .data$rank, y = .data$running_es)) +
        geom_hline(yintercept = 0, colour = "grey60", linewidth = lw) +
        geom_vline(xintercept = peak, colour = "grey60", linewidth = lw, linetype = "dashed") +
        geom_line(colour = line_col, linewidth = pt_to_linewidth(1.2)) +
        # Placed 3% of the list in from the side, so that both lines align on
        # it (a fractional hjust would shift them differently)
        annotate("text", x = if (positive) n * 0.97 else 1 + n * 0.03, y = if (positive) Inf else -Inf,
                 hjust = if (positive) 1 else 0, vjust = if (positive) 1.3 else -0.4,
                 label = stat_label, size = pt_to_size(sz$small), lineheight = 0.95, colour = "black") +
        x_scale +
        scale_y_continuous(expand = expansion(mult = c(0.06, 0.1))) +
        labs(x = NULL, y = "Running\nenrichment score") +
        stacked +
        panel_size_theme(panel.width, panel.heights[1]) +
        no_x +
        theme(plot.margin = margin(0, 0, 0, 0))

    # 2. Positions of the set's genes (ticks) above a band of the ranking
    #    metric: magenta where genes rank up, blue where they rank down, white
    #    around 0. The colours saturate at the 95th percentile of |metric|, so
    #    that a few extreme genes do not leave the rest of the band white. The
    #    band is drawn as runs of equal colour (the metric is sorted, so
    #    colours change monotonically): a few dozen rectangles, not one per
    #    gene
    band_max <- stats::quantile(abs(run$metric), 0.95, names = FALSE)
    if (!is.finite(band_max) || band_max == 0) band_max <- max(abs(run$metric), 1e-8)
    band_pal <- scales::col_numeric(c(dir_cols[["down"]], "white", dir_cols[["up"]]), domain = c(-1, 1))
    level <- round(pmax(-1, pmin(1, run$metric / band_max)) * 50) / 50
    runs <- rle(level)
    band <- data.frame(xmin = c(0, cumsum(runs$lengths)[-length(runs$lengths)]) + 0.5,
                       xmax = cumsum(runs$lengths) + 0.5, fill = band_pal(runs$values))
    hit_df <- run[run$in_set, , drop = FALSE]
    p2 <- ggplot() +
        geom_rect(aes(xmin = .data$xmin, xmax = .data$xmax, ymin = 0, ymax = 0.3, fill = I(.data$fill)),
                  data = band, colour = NA) +
        geom_segment(aes(x = .data$rank, xend = .data$rank, y = 0.36, yend = 1), data = hit_df,
                     colour = "black", linewidth = pt_to_linewidth(0.4)) +
        x_scale +
        scale_y_continuous(limits = c(0, 1), expand = expansion(0)) +
        labs(x = NULL, y = NULL) +
        stacked +
        panel_size_theme(panel.width, panel.heights[2]) +
        no_x +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(0, 0, 0, 0))

    plots <- list(p1, p2)
    # 3. Ranked metric: the same for every gene set of a comparison, so drawn
    #    as quiet context (grey area), with the rank where it crosses 0
    if (show.metric) {
        zero <- which(run$metric <= 0)[1]
        p3 <- ggplot(run, aes(x = .data$rank, y = .data$metric)) +
            geom_area(fill = "grey75", colour = NA) +
            geom_hline(yintercept = 0, colour = "grey50", linewidth = lw) +
            (if (!is.na(zero) && zero > 1) {
                geom_vline(xintercept = zero, colour = "grey50", linewidth = lw, linetype = "dashed")
            }) +
            x_scale +
            scale_y_continuous(breaks = scales::breaks_extended(3), expand = expansion(mult = 0.05)) +
            labs(x = "Rank in ordered gene list", y = metric.title) +
            stacked +
            panel_size_theme(panel.width, panel.heights[3]) +
            theme(plot.margin = margin(0, 0, 0, 0))
        plots <- c(plots, list(p3))
    } else {
        # The gene positions close the figure with the rank axis
        plots[[2]] <- p2 + labs(x = "Rank in ordered gene list") +
            theme(axis.text.x = element_text(size = sz$small * 0.9),
                  axis.ticks.x = element_line(colour = "black", linewidth = pt_to_linewidth(sz$axis)),
                  axis.line.x = element_line(colour = "black", linewidth = pt_to_linewidth(sz$axis)),
                  axis.title.x = element_text(size = sz$base))
    }

    # Title: the gene set's name; subtitle: its ID and the comparison
    if (is.null(title)) title <- info$term
    sub_parts <- c(if (show.id && info$set_id != title) info$set_id, comparison)
    transparent_bg <- is.na(bg) || identical(bg, "transparent")
    p <- patchwork::wrap_plots(plots, ncol = 1) +
        patchwork::plot_annotation(
            title = title,
            subtitle = if (length(sub_parts)) paste(sub_parts, collapse = "\n"),
            theme = theme(
                plot.title = element_text(size = sz$large, face = "bold", hjust = 0),
                plot.subtitle = element_text(size = sz$small, face = "italic", hjust = 0, lineheight = 1.05,
                                             margin = margin(t = 2, b = 4)),
                plot.background = element_rect(fill = if (transparent_bg) "transparent" else bg, colour = NA)
            )
        )
    src <- run
    attr(p, "source_data") <- src
    p
}


#' Running enrichment score of a gene set along a ranked list
#'
#' @param stats Ranking metric, named by gene
#'
#' @param set Genes of the set
#'
#' @param p Weight exponent, default is 1 (GSEA, fgsea, clusterProfiler)
#'
#' @return A data frame with one row per gene, largest metric first: `rank`,
#'         `gene`, `metric`, `running_es` and `in_set`
#'
#' @keywords internal
gsea_running_score <- function(stats,
                               set,
                               p = 1) {
    stats <- stats[!is.na(stats)]
    stats <- sort(stats, decreasing = TRUE)
    in_set <- names(stats) %in% set
    n_hit <- sum(in_set)
    if (n_hit == 0) stop("None of the gene set's genes are in the ranked list", call. = FALSE)
    if (n_hit == length(stats)) stop("Every gene of the ranked list is in the gene set", call. = FALSE)
    w <- abs(stats)^p * in_set
    if (sum(w) == 0) w <- as.numeric(in_set)
    step <- ifelse(in_set, w / sum(w), -1 / (length(stats) - n_hit))
    data.frame(rank = seq_along(stats), gene = names(stats), metric = unname(stats),
               running_es = cumsum(step), in_set = in_set)
}
