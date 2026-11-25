# Pathway enrichment results: reader and plots


#' @title Standardize pathway enrichment results
#'
#' @description Reads enrichment results from fgsea, clusterProfiler (GSEA or
#'              over-representation) or the GSEA desktop application into one
#'              table with standard column names, for
#'              \code{\link{plot_enrichment_dots}} and
#'              \code{\link{plot_enrichment_bars}}.
#'
#' @param data A data frame of results, a clusterProfiler `gseaResult` /
#'             `enrichResult` object, or a named list of either (names =
#'             comparisons)
#'
#' @param method The tool: "auto" (default) detects it from the columns; or
#'               "fgsea", "clusterProfiler" (GSEA or over-representation,
#'               told apart by their columns) or "gsea" (GSEA desktop report,
#'               e.g. read with \code{read.delim(..., check.names = FALSE)};
#'               combine its two \code{gsea_report_for_*} files with
#'               \code{rbind()} first)
#'
#' @param comparison.col For a single table holding several comparisons, the
#'                       column naming the comparison of each row. Default is
#'                       NULL
#'
#' @param clean.terms Whether to make MSigDB-style names readable (drop the
#'                    collection prefix such as "HALLMARK_" or "GOBP_",
#'                    underscores to spaces, sentence case), default is TRUE.
#'                    Other names only get a capital first letter when their
#'                    first word is all lower case (GO descriptions). The
#'                    original names are kept in `term_id`
#'
#' @details Columns read:
#' \tabular{llllll}{
#'   tool \tab term \tab nes \tab padj \tab size \tab leading_edge \cr
#'   fgsea \tab pathway \tab NES \tab padj \tab size \tab leadingEdge \cr
#'   clusterProfiler GSEA \tab Description \tab NES \tab p.adjust \tab setSize \tab core_enrichment \cr
#'   clusterProfiler ORA \tab Description \tab none \tab p.adjust \tab BgRatio (M of M/N) \tab geneID \cr
#'   GSEA desktop \tab NAME \tab NES \tab FDR q-val \tab SIZE \tab none
#' }
#' For over-representation results, `gene_ratio` is k / n from `GeneRatio`
#' ("k/n") and `count` is `Count`. For the GSEA desktop application, `padj` is
#' its FDR q-value.
#'
#' @return A data frame with columns `term`, `term_id`, `comparison` (when
#'         several), `nes`, `pvalue`, `padj`, `size`, `gene_ratio`, `count`,
#'         `leading_edge` and `type` ("gsea" or "ora"); columns a tool does not
#'         report are NA
#'
#' @examples
#' fg <- data.frame(
#'     pathway = c("HALLMARK_HYPOXIA", "HALLMARK_E2F_TARGETS"),
#'     pval = c(1e-4, 0.2), padj = c(1e-3, 0.3), NES = c(2.1, -0.8), size = c(150, 190)
#' )
#' standardize_enrichment(fg)
#' standardize_enrichment(list("KO vs WT" = fg, "Het vs WT" = fg))
#'
#' @seealso \code{\link{plot_enrichment_dots}},
#'          \code{\link{plot_enrichment_bars}}
#'
#' @export
#'
standardize_enrichment <- function(data,
                                   method = c("auto", "fgsea", "clusterProfiler", "gsea"),
                                   comparison.col = NULL,
                                   clean.terms = TRUE) {
    method <- match.arg(method)
    if (is.list(data) && !is.data.frame(data) && !isS4(data)) {
        if (length(data) == 0) stop("`data` is an empty list", call. = FALSE)
        if (is.null(names(data)) || any(!nzchar(names(data))) || anyDuplicated(names(data))) {
            stop("A list of enrichment results must have unique names, which become the comparisons",
                 call. = FALSE)
        }
        if (!is.null(comparison.col)) {
            stop("Give either a named list of tables or `comparison.col`, not both", call. = FALSE)
        }
        parts <- lapply(names(data), function(nm) {
            out <- read_enrichment_table(data[[nm]], method, where = paste0(" in '", nm, "'"))
            out$comparison <- nm
            out
        })
        res <- do.call(rbind, parts)
        res$comparison <- factor(res$comparison, levels = names(data))
    } else {
        res <- read_enrichment_table(data, method, where = "")
        if (!is.null(comparison.col)) {
            raw <- enrichment_as_data_frame(data)
            if (!(comparison.col %in% names(raw))) stop("Column '", comparison.col, "' not found", call. = FALSE)
            cmp <- raw[[comparison.col]]
            res$comparison <- if (is.factor(cmp)) droplevels(cmp) else factor(cmp, levels = unique(cmp))
        }
    }
    res$term <- if (clean.terms) clean_term_names(res$term_id) else res$term_id
    rownames(res) <- NULL
    lead <- intersect(c("term", "term_id", "comparison"), names(res))
    res[c(lead, setdiff(names(res), lead))]
}


#' A clusterProfiler result object as a data frame
#'
#' @param x A data frame, or a clusterProfiler / DOSE result object
#'
#' @return A data frame
#'
#' @keywords internal
enrichment_as_data_frame <- function(x) {
    if (isS4(x)) {
        # The as.data.frame() methods of gseaResult / enrichResult are defined
        # in DOSE, which clusterProfiler builds on
        check_installed("DOSE", reason = "to read clusterProfiler results")
        x <- as.data.frame(x)
    }
    if (!is.data.frame(x)) stop("Enrichment results must be data frames or clusterProfiler results", call. = FALSE)
    x
}


#' Read one enrichment table into standard columns
#'
#' @param data A table or result object
#'
#' @param method Tool, or "auto"
#'
#' @param where Location for messages
#'
#' @return A data frame with the standard columns
#'
#' @keywords internal
read_enrichment_table <- function(data,
                                  method,
                                  where) {
    data <- enrichment_as_data_frame(data)
    if (nrow(data) == 0) stop("The enrichment table", where, " has no rows", call. = FALSE)
    has <- function(...) all(c(...) %in% names(data))
    # GSEA desktop column names, as read with or without check.names
    gsea_q <- intersect(c("FDR q-val", "FDR.q.val"), names(data))
    gsea_p <- intersect(c("NOM p-val", "NOM.p.val"), names(data))
    if (method == "auto") {
        method <- if (has("pathway", "NES", "padj")) {
            "fgsea"
        } else if (has("Description", "p.adjust") && (has("NES") || has("GeneRatio"))) {
            "clusterProfiler"
        } else if (has("NAME", "NES") && length(gsea_q)) {
            "gsea"
        } else {
            stop("Could not detect the enrichment tool", where, " from its columns (",
                 paste(utils::head(names(data), 8), collapse = ", "), "); set `method`", call. = FALSE)
        }
    }
    col <- function(name) if (name %in% names(data)) data[[name]] else rep(NA, nrow(data))
    collapse <- function(v) {
        if (is.list(v)) vapply(v, function(g) paste(g, collapse = "/"), character(1)) else as.character(v)
    }
    n <- nrow(data)
    out <- switch(
        method,
        fgsea = data.frame(term_id = as.character(data$pathway), nes = col("NES"), pvalue = col("pval"),
                           padj = col("padj"), size = col("size"), gene_ratio = NA_real_,
                           count = NA_real_, leading_edge = collapse(col("leadingEdge")), type = "gsea"),
        clusterProfiler = if ("NES" %in% names(data)) {
            data.frame(term_id = as.character(data$Description), nes = col("NES"), pvalue = col("pvalue"),
                       padj = col("p.adjust"), size = col("setSize"), gene_ratio = NA_real_,
                       count = NA_real_, leading_edge = collapse(col("core_enrichment")), type = "gsea")
        } else {
            data.frame(term_id = as.character(data$Description), nes = NA_real_, pvalue = col("pvalue"),
                       padj = col("p.adjust"), size = ratio_part(col("BgRatio"), 1),
                       gene_ratio = ratio_part(col("GeneRatio"), 1) / ratio_part(col("GeneRatio"), 2),
                       count = col("Count"), leading_edge = collapse(col("geneID")), type = "ora")
        },
        gsea = data.frame(term_id = as.character(data$NAME), nes = col("NES"),
                          pvalue = if (length(gsea_p)) data[[gsea_p[1]]] else NA_real_,
                          padj = data[[gsea_q[1]]], size = col("SIZE"), gene_ratio = NA_real_,
                          count = NA_real_, leading_edge = NA_character_, type = "gsea")
    )
    # GSEA desktop writes "---" for missing values
    for (v in c("nes", "pvalue", "padj", "size", "gene_ratio", "count")) {
        out[[v]] <- suppressWarnings(as.numeric(out[[v]]))
    }
    out$type <- rep_len(out$type, n)
    out
}


#' Numerator or denominator of ratios written as "k/n"
#'
#' @param x Character vector of ratios
#'
#' @param i 1 for the numerator, 2 for the denominator
#'
#' @return A numeric vector
#'
#' @keywords internal
ratio_part <- function(x,
                       i) {
    parts <- strsplit(as.character(x), "/", fixed = TRUE)
    vapply(parts, function(p) if (length(p) == 2) as.numeric(p[i]) else NA_real_, numeric(1))
}


#' Readable names for MSigDB-style gene set names
#'
#' @param x Gene set names
#'
#' @return Names without the collection prefix and underscores, in sentence
#'         case, for names in upper case with underscores. Other names get a
#'         capital first letter when their first word is all lower case
#'         (e.g. GO descriptions, "extracellular matrix organization"); words
#'         with capitals inside (e.g. "mRNA") are left as they are
#'
#' @keywords internal
clean_term_names <- function(x) {
    msig <- grepl("^[A-Z0-9_]+$", x) & grepl("_", x)
    prefix <- "^(HALLMARK|GOBP|GOCC|GOMF|GO|HP|KEGG_MEDICUS|KEGG|REACTOME|WP|BIOCARTA|PID)_"
    y <- gsub("_", " ", sub(prefix, "", x[msig]))
    y <- tolower(y)
    y <- paste0(toupper(substring(y, 1, 1)), substring(y, 2))
    # Abbreviations of the Hallmark and GO names that read wrongly in lower
    # case are restored (whole words only)
    for (a in TERM_ACRONYMS) {
        y <- gsub(paste0("\\b", tolower(a), "\\b"), a, y, ignore.case = TRUE)
    }
    x[msig] <- y
    lower_first <- !msig & grepl("^[a-z]+\\b", x) & !grepl("^[a-z]+[A-Z]", x)
    x[lower_first] <- paste0(toupper(substring(x[lower_first], 1, 1)), substring(x[lower_first], 2))
    x
}


# Abbreviations kept in upper case in cleaned MSigDB names (clean_term_names())
TERM_ACRONYMS <- c("DNA", "RNA", "ATP", "MHC", "TNFA", "NFKB", "KRAS", "MYC", "E2F", "G2M", "UV",
                   "TGF", "IL2", "IL6", "JAK", "STAT3", "STAT5", "MTORC1", "PI3K", "AKT", "MTOR",
                   "P53", "WNT", "ROS")


#' Height of a row holding labels of several lines
#'
#' @param labels Labels, lines separated by "\\n"
#'
#' @param cell.height Height asked for, in inches
#'
#' @param text.size Font size in points
#'
#' @return The larger of `cell.height` and the height of the longest label
#'
#' @keywords internal
label_row_height <- function(labels,
                             cell.height,
                             text.size) {
    lines <- max(lengths(strsplit(labels, "\n", fixed = TRUE)), 1)
    if (lines == 1) return(cell.height)
    max(cell.height, (lines * 0.9 + 0.3) * text.size / 72)
}


#' Enrichment results as a standardized table
#'
#' @param enrichment Standardized results, or anything
#'                   \code{standardize_enrichment()} reads
#'
#' @param ... Arguments to \code{standardize_enrichment()}
#'
#' @keywords internal
as_enrichment <- function(enrichment,
                          ...) {
    std <- c("term", "term_id", "nes", "padj", "type")
    if (is.data.frame(enrichment) && all(std %in% names(enrichment))) return(enrichment)
    standardize_enrichment(enrichment, ...)
}


#' Wrap long labels onto several lines
#'
#' @param x Labels
#'
#' @param width Line width in characters
#'
#' @keywords internal
wrap_labels <- function(x,
                        width) {
    vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1), USE.NAMES = FALSE)
}


# Adjusted p-value colour scales -----------------------------------------
#
# `padj.scale` of plot_enrichment_dots(), plot_enrichment_grid() and
# plot_enrichment_bars() sets how adjusted p-values map to colour:
#   "log"     colour on the -log10 scale (as "-log10"), legend labelled in
#             adjusted p-values (0.01, 10^-5): the values a reader reports
#   "-log10"  colour and legend in -log10 padj
#   "linear"  colour linear in the adjusted p-value itself, the smallest
#             darkest; sets below ~1e-3 all look alike, so it suits results
#             whose adjusted p-values are of one order of magnitude
# The helpers below give the colour value of each row and the breaks, labels,
# title and direction of the colour bar.

#' Adjusted p-values of 0 replaced by the smallest non-zero one
#'
#' @param padj Adjusted p-values
#'
#' @return `padj`, without zeros
#'
#' @details Some tools report adjusted p-values of 0 (below machine
#'          precision), whose -log10 is infinite. They are drawn as the
#'          most significant non-zero value of the results (or the smallest
#'          positive double, when all are 0).
#'
#' @keywords internal
floor_padj <- function(padj) {
    zero <- !is.na(padj) & padj <= 0
    if (!any(zero)) return(padj)
    positive <- padj[!is.na(padj) & padj > 0]
    padj[zero] <- if (length(positive)) min(positive) else .Machine$double.xmin
    padj
}


#' Colour values of adjusted p-values
#'
#' @param padj Adjusted p-values
#'
#' @param padj.scale "log", "-log10" or "linear"
#'
#' @return -log10 padj (zeros floored, see \code{floor_padj()}) for "log" and
#'         "-log10", padj for "linear"
#'
#' @keywords internal
padj_colour_value <- function(padj,
                              padj.scale) {
    if (padj.scale == "linear") return(padj)
    -log10(floor_padj(padj))
}


#' Labels of adjusted p-values
#'
#' @param x Adjusted p-values (legend breaks)
#'
#' @return Plain decimals ("0.05", "0.001") when every value is at least
#'         0.001 (or 0); otherwise all in powers of ten ("10^-5",
#'         "2.5 x 10^-4"), as an expression vector, so that one legend never
#'         mixes the two
#'
#' @keywords internal
format_padj <- function(x) {
    if (all(x == 0 | x >= 1e-3)) {
        return(format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE))
    }
    labels <- lapply(x, function(v) {
        if (v == 0) return(quote(0))
        if (v == 1) return(quote(1))
        e <- floor(log10(v))
        m <- signif(v / 10^e, 2)
        # Rounding up to 10 (e.g. 9.99e-5) moves to the next power
        if (m >= 10) {
            m <- m / 10
            e <- e + 1
        }
        # The exponent as a string: drawn with the plain hyphen, as ggplot2
        # draws negative numbers on axes (plotmath would use a minus sign)
        exponent <- as.character(e)
        if (m == 1) bquote(10^.(exponent)) else bquote(.(m) %*% 10^.(exponent))
    })
    as.expression(labels)
}


#' Breaks, labels, title and direction of an adjusted p-value colour bar
#'
#' @param limits Limits of the colour scale, in the values of
#'               \code{padj_colour_value()}
#'
#' @param padj.scale "log", "-log10" or "linear"
#'
#' @return A list: `breaks` and `labels` (waiver() for "-log10"), `title`,
#'         and `reverse` (TRUE when the colour bar is flipped so that the
#'         most significant values are at its top, as on the -log10 scales)
#'
#' @keywords internal
padj_colour_guide <- function(limits,
                              padj.scale) {
    if (padj.scale == "-log10") {
        return(list(breaks = waiver(), labels = waiver(), title = quote(-log[10] ~ "padj"), reverse = FALSE))
    }
    if (padj.scale == "log") {
        # Breaks chosen among adjusted p-values, placed at their -log10:
        # over 3 or more orders of magnitude, powers of ten at round steps
        # of the exponent (10^-5, 10^-10, ...); over narrower ranges, 1-2-5
        # steps (0.01, 0.02, 0.05)
        if (diff(range(limits)) >= 3) {
            k <- scales::breaks_extended(5, Q = c(1, 2, 5))(limits)
            k <- k[k == round(k) & k >= min(limits) & k <= max(limits)]
            b <- 10^-k
        } else {
            p_range <- sort(10^-limits)
            b <- scales::breaks_log(5)(p_range)
            b <- b[b >= p_range[1] * (1 - 1e-9) & b <= p_range[2] * (1 + 1e-9)]
        }
        return(list(breaks = -log10(b), labels = format_padj(b), title = "padj", reverse = FALSE))
    }
    b <- scales::breaks_extended(4)(limits)
    b <- b[b >= limits[1] - 1e-12 & b <= limits[2] + 1e-12]
    list(breaks = b, labels = format_padj(b), title = "padj", reverse = TRUE)
}


#' @title Grid of pathway enrichment across comparisons
#'
#' @description Shows gene sets (rows) across comparisons (columns), e.g. the
#'              Hallmark sets across cell types: colour is the normalized
#'              enrichment score (GSEA) or adjusted p-value
#'              (over-representation), size the adjusted p-value or the
#'              number of genes, and non-significant results are faded. For
#'              the top gene sets of one or a few results, see
#'              \code{\link{plot_enrichment_dots}}.
#'
#' @param enrichment Results of \code{\link{standardize_enrichment}}, or
#'                   anything it reads (a table, a clusterProfiler result, or
#'                   a named list of them, one per comparison)
#'
#' @param terms Gene sets to show (names after cleaning, or the original
#'              `term_id`s): a character vector, or a named list (groups of
#'              rows). If NULL (default), the `n` most significant of each
#'              comparison, combined
#'
#' @param n Number of top gene sets per comparison when `terms` is NULL,
#'          default is 5
#'
#' @param fill What the colour shows: "nes" or "padj" (adjusted p-value, see
#'             `padj.scale`). If NULL (default), "nes" when the results have
#'             NES (GSEA), else "padj"
#'
#' @param size What the dot size shows: "padj" (larger for smaller adjusted
#'             p-values: size grows with -log10 padj), "count" or
#'             "gene_ratio". If NULL (default), "padj" when colour shows NES,
#'             else "count"
#'
#' @param padj.cutoff Adjusted p-value below which a result is significant,
#'                    default is 0.05
#'
#' @param padj.scale How adjusted p-values are coloured (with
#'                   `fill = "padj"`): "log" (default; on the -log10 scale,
#'                   the legend labelled in adjusted p-values, e.g. 0.01 and
#'                   \eqn{10^{-5}}), "-log10" (legend in -log10 padj) or
#'                   "linear" (linear in the adjusted p-value; sets below
#'                   about 0.001 then look alike, so it suits results of one
#'                   order of magnitude). The smallest adjusted p-value is
#'                   always the darkest. With `size = "padj"`, the size
#'                   legend is labelled in adjusted p-values too, unless
#'                   "-log10"
#'
#' @param significance "fade" (default; significant dots opaque with a dark
#'                     ring, the others faded), "grey" (non-significant dots
#'                     grey) or "hide"
#'
#' @param size.max Value drawn as the largest dot (larger values are drawn at
#'                 that size). If NULL (default), 10 for -log10 adjusted
#'                 p-values, else the largest value
#'
#' @param cluster.terms Whether to order gene sets by hierarchical clustering
#'                      of their colour values across comparisons (within
#'                      their groups), default is TRUE. Otherwise in the order
#'                      given, or by the first comparison's value
#'
#' @param term.replace Replacements in gene set names, applied after
#'                     cleaning: a named character vector of regular
#'                     expressions and their replacements, e.g.
#'                     `c("Extracellular matrix" = "ECM")`. Default is NULL
#'
#' @param wrap.width Gene set names are wrapped onto lines of at most this many
#'                   characters, default is 40
#'
#' @param cols Colours of the colour scale: NULL (default; "blue-magenta" for
#'             NES, "seq-purple" for adjusted p-values), a palette name
#'             or colours from low to high (from the least to the most
#'             significant for adjusted p-values)
#'
#' @param flip Whether to draw comparisons as rows and gene sets as columns,
#'             default is FALSE
#'
#' @param cell.width,cell.height Width of a column and height of a row in
#'                               inches, default is 0.18
#'
#' @param legend.position Position of the legends, default is "right"
#'
#' @param bg Background colour, default is "transparent"
#'
#' @param ... Arguments to \code{\link{standardize_enrichment}} when
#'            `enrichment` is not standardized yet (e.g. `comparison.col`)
#'
#' @return A ggplot object; its source data (\code{\link{export_source_data}})
#'         are the plotted rows of the standardized table
#'
#' @examples
#' set.seed(1)
#' sets <- paste0("HALLMARK_SET_", 1:12)
#' gsea <- function() data.frame(pathway = sets, pval = runif(12)^3, NES = rnorm(12, sd = 1.5),
#'                               size = 100)
#' res <- lapply(list(A = gsea(), B = gsea(), C = gsea()), function(d) {
#'     d$padj <- p.adjust(d$pval, "BH")
#'     d
#' })
#' plot_enrichment_grid(res, n = 4)
#'
#' @seealso \code{\link{plot_enrichment_dots}},
#'          \code{\link{plot_enrichment_bars}},
#'          \code{\link{standardize_enrichment}}
#'
#' @import ggplot2
#' @export
#'
plot_enrichment_grid <- function(
    enrichment,
    terms = NULL,
    n = 5,
    fill = NULL,
    size = NULL,
    padj.cutoff = 0.05,
    padj.scale = c("log", "-log10", "linear"),
    significance = c("fade", "grey", "hide"),
    size.max = NULL,
    cluster.terms = TRUE,
    term.replace = NULL,
    wrap.width = 40,
    cols = NULL,
    flip = FALSE,
    cell.width = 0.18,
    cell.height = 0.18,
    legend.position = c("right", "left", "top", "bottom", "none"),
    bg = "transparent",
    ...
) {
    significance <- match.arg(significance)
    padj.scale <- match.arg(padj.scale)
    legend.position <- match.arg(legend.position)
    check_de_cutoffs(padj.cutoff, 0)
    res <- as_enrichment(enrichment, ...)
    res$term <- replace_terms(res$term, term.replace)
    if (is.null(res$comparison)) res$comparison <- factor("Enrichment")
    has_nes <- any(!is.na(res$nes))
    if (is.null(fill)) fill <- if (has_nes) "nes" else "padj"
    fill <- match.arg(fill, c("nes", "padj"))
    if (fill == "nes" && !has_nes) stop("The results have no NES (over-representation); use fill = \"padj\"", call. = FALSE)
    if (is.null(size)) size <- if (fill == "nes" || all(is.na(res$count))) "padj" else "count"
    size <- match.arg(size, c("padj", "count", "gene_ratio"))
    if (size != "padj" && all(is.na(res[[size]]))) {
        stop("The results have no `", size, "` (GSEA results have no gene counts); use size = \"padj\"", call. = FALSE)
    }
    # Adjusted p-values of 0 are drawn as the smallest non-zero one
    res$neg_log10_padj <- -log10(floor_padj(res$padj))
    res$significant <- !is.na(res$padj) & res$padj < padj.cutoff

    # Gene sets: given (by cleaned name or original id), or the top `n`
    if (is.null(terms)) {
        sig <- res[res$significant, , drop = FALSE]
        if (nrow(sig) == 0) stop("No significant gene sets to pick from; give `terms`", call. = FALSE)
        sig <- sig[order(sig$padj, -abs(sig$nes)), , drop = FALSE]
        terms <- unique(unlist(lapply(split(sig$term, sig$comparison), utils::head, n)))
    }
    ft <- dot_feature_table(terms)
    by_id <- ft$feature %in% res$term_id & !(ft$feature %in% res$term)
    ft$feature[by_id] <- res$term[match(ft$feature[by_id], res$term_id)]
    missing <- setdiff(ft$feature, res$term)
    if (length(missing)) {
        warning(length(missing), " gene set(s) not in the results are left out: ",
                paste(missing, collapse = ", "), call. = FALSE)
        ft <- ft[!(ft$feature %in% missing), , drop = FALSE]
    }
    if (nrow(ft) == 0) stop("None of `terms` are in the results", call. = FALSE)
    ft$id <- if (is.null(ft$feature_group)) ft$feature else paste(ft$feature_group, ft$feature, sep = "\r")
    labels <- stats::setNames(wrap_labels(ft$feature, wrap.width), ft$id)
    # Rows grow to hold names wrapped onto several lines
    if (!flip) cell.height <- label_row_height(labels, cell.height, theme_sizes()$small)

    # One row per gene set and comparison
    cmp <- levels(res$comparison)
    idx <- expand.grid(j = seq_len(nrow(ft)), i = seq_along(cmp))
    df <- data.frame(id = ft$id[idx$j], term = ft$feature[idx$j], comparison = cmp[idx$i],
                     stringsAsFactors = FALSE)
    if (!is.null(ft$feature_group)) df$term_group <- ft$feature_group[idx$j]
    hit <- match(paste(df$comparison, df$term, sep = "\r"), paste(res$comparison, res$term, sep = "\r"))
    for (v in c("nes", "pvalue", "padj", "neg_log10_padj", "size", "gene_ratio", "count", "significant")) {
        df[[v]] <- res[[v]][hit]
    }
    df$significant <- df$significant %in% TRUE
    df$comparison <- factor(df$comparison, levels = cmp)

    # Colour: NES, or the adjusted p-value on the scale of `padj.scale`
    fill_v <- if (fill == "nes") df$nes else padj_colour_value(df$padj, padj.scale)
    size_v <- switch(size, padj = df$neg_log10_padj, count = df$count, gene_ratio = df$gene_ratio)
    drawn <- !is.na(df$padj) & !is.na(size_v)
    if (significance == "hide") drawn <- drawn & df$significant
    if (significance == "grey") fill_v[!df$significant] <- NA

    # Order of gene sets: clustered on the colour values (adjusted p-values
    # on the -log10 scale, whatever `padj.scale`), or as given
    term_lv <- unique(ft$id)
    if (cluster.terms && length(term_lv) > 2) {
        cluster_v <- if (fill == "nes" || padj.scale != "linear") fill_v else df$neg_log10_padj
        if (significance == "grey") cluster_v[!df$significant] <- NA
        m <- tapply(cluster_v, list(df$id, as.character(df$comparison)), function(v) v[1])
        term_lv <- cluster_order(m[term_lv, cmp, drop = FALSE], ft$feature_group)
    }

    if (is.null(size.max)) {
        size.max <- if (size == "padj") 10 else max(size_v[drawn], na.rm = TRUE)
    }
    breaks <- scales::breaks_extended(4)(c(0, size.max))
    breaks <- breaks[breaks > 0 & breaks <= size.max]
    capped <- any(size_v[drawn] > size.max, na.rm = TRUE)
    if (capped && utils::tail(breaks, 1) != size.max) breaks <- c(breaks, size.max)
    # Sizes grow with -log10 padj whatever `padj.scale`; unless "-log10",
    # their legend is labelled in adjusted p-values, as the colour bar
    padj_size <- size == "padj" && padj.scale != "-log10"
    blabels <- if (padj_size) format_padj(10^-breaks) else format(breaks, trim = TRUE)
    if (capped) {
        # The largest dot stands for this value and beyond: at least this
        # -log10 padj, at most this adjusted p-value
        sign <- if (padj_size) "\u2264" else "\u2265"
        k <- length(blabels)
        if (is.expression(blabels)) {
            blabels[[k]] <- bquote(.(sign) * .(blabels[[k]]))
        } else {
            blabels[k] <- paste0(sign, blabels[k])
        }
    }
    fv <- fill_v[drawn & !is.na(fill_v)]
    limits <- if (fill == "nes") {
        # Centred on 0, each side to its own extreme (diverging_limits())
        diverging_limits(fv, center = 0)
    } else if (padj.scale == "linear") {
        c(0, max(c(fv, padj.cutoff)))
    } else {
        c(0, max(c(fv, -log10(padj.cutoff))))
    }
    fill_cols <- resolve_continuous_palette(cols, diverging = fill == "nes", default = "seq-purple")
    fill_guide <- if (fill == "nes") {
        list(breaks = waiver(), labels = waiver(), title = "NES", reverse = FALSE)
    } else {
        padj_colour_guide(limits, padj.scale)
    }
    # Linear adjusted p-values: the smallest (most significant) darkest
    if (fill == "padj" && padj.scale == "linear") fill_cols <- rev(fill_cols)

    term_f <- factor(df$id, levels = term_lv)
    cmp_f <- df$comparison
    grp <- if (!is.null(ft$feature_group)) factor(df$term_group, levels = levels(ft$feature_group))
    eng <- data.frame(size = ifelse(drawn, size_v, NA), fill = fill_v)
    if (!flip) {
        eng$x <- cmp_f
        eng$y <- term_f
        eng$y_group <- grp
    } else {
        eng$x <- term_f
        eng$y <- cmp_f
        eng$x_group <- grp
    }
    if (significance == "fade") eng$outline <- df$significant
    cmp_labels <- stats::setNames(cmp, cmp)
    p <- dot_grid_plot(
        eng,
        x.labels = if (!flip) cmp_labels else labels,
        y.labels = if (!flip) labels else cmp_labels,
        size.limits = c(0, size.max),
        size.breaks = breaks,
        size.labels = blabels,
        size.title = switch(
            size,
            padj = if (padj_size) "padj" else quote(-log[10] ~ "padj"),
            count = "Genes",
            gene_ratio = "Gene ratio"
        ),
        fill.colours = fill_cols,
        fill.limits = limits,
        fill.center = if (fill == "nes") 0 else NULL,
        fill.title = fill_guide$title,
        fill.breaks = fill_guide$breaks,
        fill.labels = fill_guide$labels,
        fill.reverse = fill_guide$reverse,
        na.value = "grey80",
        outline.labels = significance_labels(padj.cutoff),
        outline.title = "Significance",
        cell.width = cell.width,
        cell.height = cell.height,
        legend.position = legend.position,
        bg = bg
    )
    keep <- intersect(c("term_group", "term", "comparison", "nes", "pvalue", "padj", "size", "gene_ratio",
                        "count", "significant"), names(df))
    attr(p, "source_data") <- df[, keep]
    p
}


#' Replace patterns in gene set names
#'
#' @param x Gene set names
#'
#' @param replace A named character vector: regular expressions (names) and
#'                their replacements, applied in order; or NULL
#'
#' @return The names with the replacements made
#'
#' @keywords internal
replace_terms <- function(x,
                          replace) {
    if (is.null(replace)) return(x)
    if (!(is.character(replace) && !is.null(names(replace)) && all(nzchar(names(replace))))) {
        stop("`term.replace` must be a named character vector: c(pattern = \"replacement\")", call. = FALSE)
    }
    for (i in seq_along(replace)) x <- gsub(names(replace)[i], replace[[i]], x, perl = TRUE)
    x
}


#' @title Dot plot of the top gene sets of enrichment results
#'
#' @description Shows the most significant gene sets of one enrichment result,
#'              or of several side by side (e.g. the gene sets enriched in up-
#'              and in down-regulated genes): the gene ratio (or NES) on the
#'              x-axis, the number of genes as dot size and the adjusted
#'              p-value as colour (darker for smaller). One result is drawn in
#'              reversed magma; several are panels named by strips, each
#'              drawn in its own single-hue palette, so that panels read apart
#'              at a glance; panels are stacked (default), placed side by side,
#'              or in a grid. Panels have a fixed size, so the figure combines
#'              cleanly with others (e.g. two timepoints with
#'              \code{patchwork}); save it with \code{\link{save_plot}}.
#'
#' @param enrichment One enrichment result (a table or a clusterProfiler
#'                   result, see \code{\link{standardize_enrichment}}), or a
#'                   named list of them, one panel each, e.g.
#'                   `list("Up-regulated" = ora_up, "Down-regulated" = ora_down)`.
#'                   Standardized results with a `comparison` column give one
#'                   panel per comparison
#'
#' @param n Number of gene sets per panel, the most significant first, default
#'          is 10
#'
#' @param terms Gene sets to show instead of the top `n` (names after
#'              cleaning, or the original `term_id`s); each panel shows those
#'              it has. Default is NULL
#'
#' @param x What the x-axis shows: "gene_ratio" (default; over-representation:
#'          genes of the set in the list / genes in the list, `k/n`; GSEA:
#'          leading-edge genes / genes in the set), "nes" (GSEA), "count" or
#'          "padj" (-log10 adjusted p-value)
#'
#' @param order.by Order of the gene sets in a panel, from the top: "x"
#'                 (default; largest x first) or "padj" (most significant
#'                 first)
#'
#' @param padj.cutoff Adjusted p-value below which a gene set is significant,
#'                    default is 0.05
#'
#' @param padj.scale How adjusted p-values are coloured: "log" (default;
#'                   on the -log10 scale, the legend labelled in adjusted
#'                   p-values, e.g. 0.01 and \eqn{10^{-5}}), "-log10" (legend
#'                   in -log10 padj) or "linear" (linear in the adjusted
#'                   p-value; sets below about 0.001 then look alike, so it
#'                   suits results of one order of magnitude). The smallest
#'                   adjusted p-value is always the darkest
#'
#' @param show.ns Whether non-significant gene sets may fill a panel that has
#'                fewer than `n` significant ones, default is FALSE
#'
#' @param term.replace Replacements in gene set names, applied after
#'                     cleaning: a named character vector of regular
#'                     expressions and their replacements, e.g.
#'                     `c("Extracellular matrix" = "ECM",
#'                     "transforming growth factor beta" = "TGF-beta")`;
#'                     Greek letters as Unicode escapes, e.g. "TGF\\u03b2"
#'                     for TGF-beta with the Greek letter.
#'                     Default is NULL
#'
#' @param wrap.width Gene set names are wrapped onto lines of at most this
#'                   many characters, default is 45
#'
#' @param cols Colours of the adjusted p-value, from the least to the most
#'             significant. NULL (default) gives one result reversed magma
#'             (light orange to near-black purple, trimmed at both ends), and
#'             several results each their own single-hue palette, in the
#'             order "seq-purple", "seq-blue", "seq-green", "seq-orange",
#'             "seq-red", "seq-teal", "seq-magenta"; one palette name (or
#'             colours) for every panel; or palette names, one per panel
#'             (named by panel, or in panel order), or a list of colour
#'             vectors
#'
#' @param shared.scale Whether every panel uses the same colour and size
#'                     limits (and x-axis range), default is TRUE, so that the
#'                     same lightness and dot size mean the same value in every
#'                     panel. Set FALSE to give each panel its own
#'
#' @param ncol,nrow Layout of the panels. If both are NULL (default), panels
#'                  are stacked in one column; `nrow = 1` places them side by
#'                  side
#'
#' @param title Title of the figure. Default is NULL
#'
#' @param fill.title,size.title Titles of the colour and size legends. If NULL
#'                              (default), "padj" (-log10 padj with
#'                              `padj.scale = "-log10"`) and "Genes"
#'
#' @param cell.height Height of a gene set's row in inches, default is 0.2
#'                    (rows of names wrapped onto several lines are taller)
#'
#' @param panel.width Width of a panel in inches, default is 2
#'
#' @param legend.position Position of the legends of each panel: "right"
#'                        (default), "bottom" or "none"
#'
#' @param legend.layout Arrangement of the legends on the right:
#'                      "columns" (default; colour bars in one column and
#'                      size legends in the next, each in panel order; a
#'                      single panel's two legends side by side) or "stacked"
#'                      (one column, each panel's colour bar and size legend
#'                      in turn)
#'
#' @param bg Background colour, default is "transparent"
#'
#' @param ... Arguments to \code{\link{standardize_enrichment}} when
#'            `enrichment` is not standardized yet (e.g. `clean.terms`)
#'
#' @details Legends: each panel has its own colour bar (in its palette) when
#'          panels are drawn in different palettes; the size legend is drawn
#'          once when the limits are shared. With legends on the right
#'          (default), the legends of all panels are gathered to the right of
#'          the panels: colour bars in one column and size legends in the
#'          next (`legend.layout = "columns"`), or all in one column. Panels
#'          side by side (or in a grid) have the same height, so that their
#'          rows line up: shorter panels get empty rows at the bottom. Long
#'          names of stacked panels are wrapped to fit the height of their
#'          strip; all strips have the same thickness, with their names
#'          centred.
#'
#'          Gene ratio is a descriptive measure: a set of 10 genes with 5 in
#'          the list has a higher ratio than a set of 500 with 50, but is
#'          supported by fewer genes; read it with the dot size (count) and
#'          colour (significance).
#'
#' @return A patchwork of one panel per result (a ggplot for a single
#'         result); its source data (\code{\link{export_source_data}}) are the
#'         plotted rows of the standardized table, with their `panel`
#'
#' @examples
#' set.seed(1)
#' ora <- function(k) {
#'     d <- data.frame(Description = paste("gene set", 1:15), Count = k,
#'                     GeneRatio = paste0(k, "/200"), BgRatio = "100/20000",
#'                     pvalue = runif(15)^6)
#'     d$p.adjust <- p.adjust(d$pvalue, "BH")
#'     d
#' }
#' plot_enrichment_dots(ora(sample(5:40, 15)))
#' plot_enrichment_dots(list("Up-regulated" = ora(sample(5:40, 15)),
#'                           "Down-regulated" = ora(sample(5:40, 15))),
#'                      n = 6, title = "GO biological process")
#'
#' @seealso \code{\link{plot_enrichment_grid}},
#'          \code{\link{plot_enrichment_bars}},
#'          \code{\link{plot_gsea_rank}}, \code{\link{standardize_enrichment}}
#'
#' @import ggplot2
#' @export
#'
plot_enrichment_dots <- function(
    enrichment,
    n = 10,
    terms = NULL,
    x = c("gene_ratio", "nes", "count", "padj"),
    order.by = c("x", "padj"),
    padj.cutoff = 0.05,
    padj.scale = c("log", "-log10", "linear"),
    show.ns = FALSE,
    term.replace = NULL,
    wrap.width = 45,
    cols = NULL,
    shared.scale = TRUE,
    ncol = NULL,
    nrow = NULL,
    title = NULL,
    fill.title = NULL,
    size.title = NULL,
    cell.height = 0.2,
    panel.width = 2,
    legend.position = c("right", "bottom", "none"),
    legend.layout = c("columns", "stacked"),
    bg = "transparent",
    ...
) {
    x <- match.arg(x)
    order.by <- match.arg(order.by)
    padj.scale <- match.arg(padj.scale)
    legend.position <- match.arg(legend.position)
    legend.layout <- match.arg(legend.layout)
    check_de_cutoffs(padj.cutoff, 0)
    res <- as_enrichment(enrichment, ...)
    res$term <- replace_terms(res$term, term.replace)
    # A single result is one panel without a strip
    has_strip <- !is.null(res$comparison)
    if (!has_strip) res$comparison <- factor("Enrichment")
    res$comparison <- droplevels(res$comparison)
    # One result (a table, or a list of one) is drawn in reversed magma;
    # counted before panels without gene sets are left out, so that the
    # palettes do not change with the results
    single_result <- nlevels(res$comparison) == 1
    res <- res[!is.na(res$padj), , drop = FALSE]
    # Adjusted p-values of 0 are drawn as the smallest non-zero one
    res$neg_log10_padj <- -log10(floor_padj(res$padj))
    # Colour value: -log10 padj ("log", "-log10") or padj ("linear")
    res$padj_colour <- padj_colour_value(res$padj, padj.scale)
    res$significant <- res$padj < padj.cutoff

    # GSEA: count and gene ratio from the leading edge
    gsea <- res$type == "gsea"
    if (any(gsea)) {
        le <- ifelse(is.na(res$leading_edge[gsea]) | res$leading_edge[gsea] == "", NA_character_,
                     res$leading_edge[gsea])
        res$count[gsea] <- ifelse(is.na(le), NA, lengths(strsplit(le, "/", fixed = TRUE)))
        res$gene_ratio[gsea] <- res$count[gsea] / res$size[gsea]
    }
    res$x <- switch(x, gene_ratio = res$gene_ratio, nes = res$nes, count = res$count, padj = res$neg_log10_padj)
    if (all(is.na(res$x))) {
        stop("The results have no ", switch(x, gene_ratio = "gene ratio (GSEA results need leading-edge genes)",
                                             nes = "NES (over-representation results)", count = "gene counts",
                                             padj = "adjusted p-values"),
             "; choose another `x`", call. = FALSE)
    }
    has_count <- !all(is.na(res$count))

    # Gene sets of each panel: given, or the top `n` by adjusted p-value
    panels <- levels(res$comparison)
    picked <- lapply(panels, function(pn) {
        d <- res[res$comparison == pn & !is.na(res$x), , drop = FALSE]
        if (!is.null(terms)) {
            d <- d[d$term %in% terms | d$term_id %in% terms, , drop = FALSE]
        } else {
            if (!show.ns) d <- d[d$significant, , drop = FALSE]
            d <- utils::head(d[order(d$padj, -abs(d$x)), , drop = FALSE], n)
        }
        if (nrow(d) == 0) {
            message("No ", if (is.null(terms)) "significant " else "selected ", "gene sets in '", pn,
                    "'; its panel is left out", if (is.null(terms) && !show.ns) " (see `show.ns`)")
            return(NULL)
        }
        # Top to bottom: largest x (or smallest padj) first
        d <- if (order.by == "x") d[order(-d$x, d$padj), , drop = FALSE] else d[order(d$padj, -d$x), , drop = FALSE]
        d
    })
    names(picked) <- panels
    picked <- Filter(Negate(is.null), picked)
    if (length(picked) == 0) stop("No gene sets to draw", call. = FALSE)

    # Colours: one palette per panel
    pals <- enrichment_panel_palettes(cols, names(picked), single = single_result)

    # Shared limits: colour, size and x range over all panels
    all_d <- do.call(rbind, picked)
    range_of <- function(v) {
        r <- range(v, na.rm = TRUE)
        if (diff(r) == 0) r + c(-0.5, 0.5) * max(abs(r[1]) * 0.1, 0.1) else r
    }
    sz <- theme_sizes()
    max_size <- 0.85 * min(cell.height, 0.25) * 72 / (0.75 * ggplot2::.pt)
    same_palette <- length(unique(vapply(pals, paste, "", collapse = ""))) == 1
    if (is.null(fill.title)) fill.title <- if (padj.scale == "-log10") quote(-log[10] ~ "padj") else "padj"
    size.title <- if (is.null(size.title)) "Genes" else size.title
    x_title <- switch(x, gene_ratio = "Gene ratio", nes = "Normalized enrichment score", count = "Genes",
                      padj = quote(-log[10] ~ "padj"))

    # Panels stacked in one column (default) get vertical strips on the right
    vertical <- is.null(ncol) && is.null(nrow) || identical(ncol, 1) || identical(ncol, 1L)
    # Panels side by side (or in a grid) get the same height, so that their
    # rows and names line up: as many rows as the longest panel (empty rows
    # at the bottom of shorter ones), each as tall as the tallest
    row_hs <- vapply(picked, function(d) label_row_height(wrap_labels(d$term, wrap.width), cell.height, sz$small),
                     numeric(1))
    n_rows_max <- max(vapply(picked, nrow, integer(1)))

    # Strip names: those of stacked panels are wrapped to the height of their
    # panel, so that a long name on a short panel stays within its strip. All
    # strips get the thickness of the one with the most lines, the names
    # centred in it
    strip_names <- vapply(seq_along(picked), function(i) {
        pn <- names(picked)[i]
        if (!vertical) return(pn)
        n_rows <- nrow(picked[[i]])
        fit <- floor(grid_panel_size(n_rows, row_hs[[i]]) / (0.6 * sz$small / 72))
        paste(strwrap(pn, width = max(fit, 6)), collapse = "\n")
    }, character(1))
    # Thickness of each name as drawn (bold, 0.9 line height), measured on a
    # PDF device; each strip is padded by half its difference to the thickest
    grDevices::pdf(NULL)
    name_h <- tryCatch(vapply(strip_names, function(nm) {
        g <- grid::textGrob(nm, gp = grid::gpar(fontsize = sz$small, fontface = "bold", lineheight = 0.9))
        grid::convertHeight(grid::grobHeight(g), "pt", valueOnly = TRUE)
    }, numeric(1)), finally = grDevices::dev.off())
    strip_pad <- unname(max(name_h) - name_h) / 2
    # Legends shown by each panel: colour bars first, size second. With
    # shared limits, a colour bar is repeated only when the palettes differ
    # and the size legend is drawn once (first panel)
    show_fill <- !shared.scale | !same_palette | seq_along(picked) == 1
    show_size <- has_count & (!shared.scale | seq_along(picked) == 1)

    plots <- lapply(seq_along(picked), function(i) {
        d <- picked[[i]]
        pn <- names(picked)[i]
        lim_d <- if (shared.scale) all_d else d
        d$label <- wrap_labels(d$term, wrap.width)
        labels <- stats::setNames(d$label, d$term_id)
        row_h <- if (vertical) row_hs[[i]] else max(row_hs)
        n_rows <- if (vertical) nrow(d) else n_rows_max
        # Discrete y levels from the bottom: padding rows, then the gene sets
        # (the first at the top)
        pad <- if (n_rows > nrow(d)) paste0(".pad", seq_len(n_rows - nrow(d))) else character(0)
        y_levels <- c(pad, rev(unique(d$term_id)))
        d$id <- factor(d$term_id, levels = y_levels)
        d$panel <- factor(strip_names[i], levels = strip_names[i])
        pal <- pals[[pn]]
        # Colour scale of the adjusted p-values: limits, breaks and labels
        # (see padj_colour_guide()); linear adjusted p-values run the palette
        # backwards, so that the smallest is the darkest
        colour_lim <- range_of(lim_d$padj_colour)
        if (padj.scale == "linear") colour_lim[1] <- max(colour_lim[1], 0)
        colour_guide <- padj_colour_guide(colour_lim, padj.scale)
        scale_pal <- if (padj.scale == "linear") rev(pal) else pal
        fill_guide <- if (show_fill[i]) {
            guide_colourbar(order = 1, reverse = colour_guide$reverse, theme = theme(
                legend.key.height = unit(min(0.9, max(0.5, nrow(d) * row_h * 0.8)), "in"),
                legend.key.width = unit(0.12, "in"),
                legend.frame = element_rect(colour = "black", linewidth = 0.3),
                legend.ticks = element_line(colour = "black", linewidth = 0.3),
                legend.ticks.length = unit(0.03, "in")
            ))
        } else {
            "none"
        }
        # Solid dots without an outline (a ring would darken the pale end of
        # the palette); without counts (GSEA desktop reports) all of one size
        dots <- if (has_count) {
            geom_point(aes(colour = .data$padj_colour, size = .data$count), shape = 16)
        } else {
            geom_point(aes(colour = .data$padj_colour), shape = 16, size = max_size * 0.6)
        }
        p <- ggplot(d, aes(x = .data$x, y = .data$id)) +
            dots +
            scale_colour_gradientn(
                colours = scale_pal,
                limits = colour_lim,
                breaks = colour_guide$breaks,
                labels = colour_guide$labels,
                oob = scales::oob_squish,
                name = bold_title(fill.title),
                guide = fill_guide
            ) +
            # Padding rows have no tick, grid line or name
            scale_y_discrete(limits = y_levels, breaks = setdiff(y_levels, pad), labels = labels,
                             expand = grid_expansion()) +
            scale_x_continuous(limits = if (shared.scale) range(all_d$x) else NULL,
                               expand = expansion(mult = 0.12)) +
            labs(x = x_title, y = NULL) +
            theme_bioToolkit("axis", bg = bg) +
            panel_size_theme(panel.width, grid_panel_size(n_rows, row_h)) +
            theme(
                axis.line = element_blank(),
                panel.border = element_rect(colour = "black", fill = NA, linewidth = pt_to_linewidth(sz$line)),
                panel.grid.major = element_line(colour = "grey90", linewidth = pt_to_linewidth(sz$line)),
                axis.ticks = element_line(colour = "black", linewidth = pt_to_linewidth(sz$line)),
                axis.text.y = element_text(lineheight = 0.85),
                legend.position = legend.position,
                legend.title = element_text(size = sz$small, margin = margin(b = 0.15, unit = "cm")),
                legend.title.position = "top",
                legend.spacing.y = unit(0.15, "in"),
                legend.justification = c(0, 1)
            )
        if (has_count) {
            p <- p + scale_size(range = c(max_size * 0.3, max_size), limits = range_of(lim_d$count),
                                name = bold_title(size.title),
                                guide = if (show_size[i]) {
                                    guide_legend(order = 2, override.aes = list(colour = "grey45"))
                                } else {
                                    "none"
                                })
        }
        if (has_strip) {
            # The strip is tinted with (half of) the light end of the panel's
            # palette and named in its dark end, tying the panel to its colour
            # bar
            # Padded across the strip to the thickness of the thickest one
            pad <- 3 + strip_pad[i]
            strip_theme <- element_text(
                size = sz$small,
                face = "bold",
                colour = pal[length(pal)],
                lineheight = 0.9,
                margin = if (vertical) margin(3, pad, 3, pad) else margin(pad, 3, pad, 3)
            )
            p <- p +
                (if (vertical) facet_grid(rows = vars(.data$panel)) else facet_grid(cols = vars(.data$panel))) +
                theme(
                    strip.background = element_rect(fill = grDevices::colorRampPalette(c("white", pal[1]))(3)[2],
                                                    colour = "black",
                                                    linewidth = pt_to_linewidth(sz$line)),
                    strip.text.y.right = strip_theme,
                    strip.text.x.top = strip_theme,
                    strip.clip = "off"
                )
        }
        p
    })

    src <- do.call(rbind, lapply(names(picked), function(pn) {
        d <- picked[[pn]]
        data.frame(panel = pn, term = d$term, term_id = d$term_id, x = d$x, gene_ratio = d$gene_ratio,
                   count = d$count, nes = d$nes, pvalue = d$pvalue, padj = d$padj, size = d$size,
                   significant = d$significant)
    }))
    names(src)[names(src) == "x"] <- paste0("x_", x)
    if (length(plots) == 1) {
        p <- plots[[1]]
        # Colour bar and size legend side by side
        if (legend.layout == "columns") p <- p + theme(legend.box = "horizontal", legend.box.just = "top")
        if (!is.null(title)) p <- p + labs(title = title) + theme(plot.title.position = "plot")
        attr(p, "source_data") <- src
        return(p)
    }

    # Panels stacked by default; with their x-axes kept, each panel reads on
    # its own when cut out of the figure
    if (is.null(ncol) && is.null(nrow)) ncol <- 1
    transparent_bg <- is.na(bg) || identical(bg, "transparent")
    # Legends are collected into one column on the right (each panel's colour
    # bar in its palette, the size legend once): kept per panel, the legends
    # of a panel shorter than them would spill over the next panel or out of
    # the figure
    guides <- if (legend.position == "right") "collect" else "keep"
    # Stacked panels with a shared x range keep only the bottom x-axis, as
    # facets would (patchwork drops axes identical to their neighbour's)
    axes <- if (shared.scale && identical(as.numeric(ncol), 1)) "collect_x" else "keep"
    annotation <- patchwork::plot_annotation(title = title, theme = theme(
        plot.title = element_text(size = sz$large, face = "bold", hjust = 0),
        plot.background = element_rect(fill = if (transparent_bg) "transparent" else bg, colour = NA)
    ))
    if (legend.position == "right" && legend.layout == "columns") {
        # Legends in two columns to the right of the panels: the colour bars,
        # then the size legends, each in panel order. Each legend is taken
        # from its panel (the other legend hidden) when the figure is drawn
        fill_col <- lapply(which(show_fill), function(i) plots[[i]] + guides(size = "none"))
        size_col <- lapply(which(show_size), function(i) plots[[i]] + guides(colour = "none"))
        block <- legend_block_grob(list(fill_col, size_col))
        legend_fig <- legend_block_figure(block)
        panels <- lapply(plots, function(p) p + theme(legend.position = "none"))
        # One flat layout (a nested patchwork of fixed-size panels cannot be
        # printed, a patchwork limitation): the panels in their grid, the
        # legends in a column to the right spanning all rows, top aligned
        grid_ncol <- if (!is.null(ncol)) ncol else ceiling(length(panels) / nrow)
        grid_nrow <- ceiling(length(panels) / grid_ncol)
        k <- seq_along(panels) - 1
        design <- do.call(c, lapply(seq_along(panels), function(i) {
            patchwork::area(t = k[i] %/% grid_ncol + 1, l = k[i] %% grid_ncol + 1)
        }))
        # Legends taller than the rows of panels get a blank row below them
        grDevices::pdf(NULL)
        span <- tryCatch({
            gt <- patchwork::patchworkGrob(patchwork::wrap_plots(panels, ncol = grid_ncol, axes = axes))
            pos <- gt$layout[grepl("^panel", gt$layout$name), , drop = FALSE]
            sum(grid::convertHeight(gt$heights, "in", valueOnly = TRUE)[min(pos$t):max(pos$b)])
        }, finally = grDevices::dev.off())
        deficit <- attr(legend_fig, "block_height") - span
        bottom <- grid_nrow
        if (deficit > 0) {
            bottom <- grid_nrow + 1
            spacer <- ggplot() +
                theme_void() +
                panel_size_theme(0.01, deficit) +
                theme(plot.margin = margin(0, 0, 0, 0))
            panels <- c(panels, list(spacer))
            design <- c(design, patchwork::area(t = bottom, l = 1))
        }
        design <- c(design, patchwork::area(t = 1, l = grid_ncol + 1, b = bottom))
        combined <- patchwork::wrap_plots(c(panels, list(legend_fig)), design = design, axes = axes) +
            annotation
    } else {
        combined <- patchwork::wrap_plots(plots, ncol = ncol, nrow = nrow, guides = guides, axes = axes) +
            annotation
    }
    attr(combined, "source_data") <- src
    combined
}


#' Palettes of the panels of an enrichment dot plot
#'
#' @param cols NULL; one palette name or colour vector for all panels;
#'             palette names, one per panel; or a list of colour vectors
#'
#' @param panels Panel names
#'
#' @param single Whether the plot has one result; with `cols = NULL`, it is
#'               drawn in reversed magma (\code{magma_r_trimmed()}). Default
#'               is FALSE
#'
#' @return A list of colour vectors, named by panel
#'
#' @keywords internal
enrichment_panel_palettes <- function(cols,
                                      panels,
                                      single = FALSE) {
    if (is.null(cols) && single) return(stats::setNames(list(magma_r_trimmed()), panels))
    defaults <- c("seq-purple", "seq-blue", "seq-green", "seq-orange", "seq-red", "seq-teal", "seq-magenta")
    per_panel <- if (is.null(cols)) {
        as.list(rep_len(defaults, length(panels)))
    } else if (is.list(cols)) {
        cols
    } else if (length(cols) > 1 && all(vapply(cols, function(cl) !is.null(find_palette(cl)), logical(1)))) {
        as.list(cols)
    } else {
        rep(list(cols), length(panels))
    }
    if (!is.null(names(per_panel))) {
        missing <- setdiff(panels, names(per_panel))
        if (length(missing)) stop("No colours given for panel(s): ", paste(missing, collapse = ", "), call. = FALSE)
        per_panel <- per_panel[panels]
    } else {
        if (length(per_panel) < length(panels)) {
            stop("`cols` has ", length(per_panel), " palettes for ", length(panels), " panels", call. = FALSE)
        }
        per_panel <- per_panel[seq_along(panels)]
    }
    stats::setNames(lapply(per_panel, resolve_continuous_palette, diverging = FALSE), panels)
}


#' @title Diverging bar chart of normalized enrichment scores
#'
#' @description Shows the most significant gene sets of one GSEA comparison as
#'              bars of their normalized enrichment score (NES), positive to
#'              the right and negative to the left of 0, ordered by NES. Each
#'              gene set's name is written on the opposite side of 0 from its
#'              bar, so that names never cover bars and no axis labels are
#'              needed.
#'
#' @param enrichment Results of \code{\link{standardize_enrichment}}, or
#'                   anything it reads, with NES (GSEA results)
#'
#' @param comparison The comparison to draw, when the results hold several.
#'                   If NULL (default), the first
#'
#' @param n.up,n.down Number of gene sets with positive / negative NES to
#'                    show, the most significant first, default is 10 each
#'
#' @param fill Bar colours: "direction" (default; one colour for positive and
#'             one for negative NES) or "padj" (adjusted p-value, see
#'             `padj.scale`)
#'
#' @param padj.cutoff Adjusted p-value below which a gene set is significant,
#'                    default is 0.05
#'
#' @param padj.scale How adjusted p-values are coloured (with
#'                   `fill = "padj"`): "log" (default; on the -log10 scale,
#'                   the legend labelled in adjusted p-values, e.g. 0.01 and
#'                   \eqn{10^{-5}}), "-log10" (legend in -log10 padj) or
#'                   "linear" (linear in the adjusted p-value). The smallest
#'                   adjusted p-value is always the darkest
#'
#' @param show.ns Whether non-significant gene sets may be shown (faded, and
#'                explained in a "Significance" legend), when there are
#'                fewer than `n.up` / `n.down` significant ones, default is
#'                FALSE
#'
#' @param cols With `fill = "direction"`, the colours of positive and
#'             negative NES: a vector named `up` and / or `down`, or two
#'             colours (up, down). NULL (default) takes the ends of the
#'             "blue-magenta" palette. With `fill = "padj"`, a sequential
#'             palette name or colours from the least to the most
#'             significant, default "seq-purple"
#'
#' @param term.replace Replacements in gene set names, applied after
#'                     cleaning: a named character vector of regular
#'                     expressions and their replacements. Default is NULL
#'
#' @param wrap.width Names are wrapped onto lines of at most this many
#'                   characters. If NULL (default), as wide as half the panel
#'                   allows
#'
#' @param cell.height Height of a bar's row in inches, default is 0.2
#'
#' @param panel.width Width of the panel in inches, default is 4
#'
#' @param legend.position Position of the legends, default is "right". With
#'                        `fill = "direction"` the colours need no legend (the
#'                        side of 0 tells the direction); the significance
#'                        legend is drawn when non-significant bars are shown
#'
#' @param bg Background colour, default is "transparent"
#'
#' @param ... Arguments to \code{\link{standardize_enrichment}}
#'
#' @return A ggplot object; its data are the plotted gene sets
#'
#' @examples
#' set.seed(1)
#' fg <- data.frame(pathway = paste0("HALLMARK_SET_", 1:30), pval = runif(30)^4,
#'                  NES = rnorm(30, sd = 1.5), size = 100)
#' fg$padj <- p.adjust(fg$pval, "BH")
#' plot_enrichment_bars(fg, n.up = 5, n.down = 5)
#'
#' @seealso \code{\link{plot_enrichment_dots}},
#'          \code{\link{standardize_enrichment}}
#'
#' @import ggplot2
#' @export
#'
plot_enrichment_bars <- function(
    enrichment,
    comparison = NULL,
    n.up = 10,
    n.down = 10,
    fill = c("direction", "padj"),
    padj.cutoff = 0.05,
    padj.scale = c("log", "-log10", "linear"),
    show.ns = FALSE,
    cols = NULL,
    term.replace = NULL,
    wrap.width = NULL,
    cell.height = 0.2,
    panel.width = 4,
    legend.position = "right",
    bg = "transparent",
    ...
) {
    fill <- match.arg(fill)
    padj.scale <- match.arg(padj.scale)
    check_de_cutoffs(padj.cutoff, 0)
    res <- as_enrichment(enrichment, ...)
    res$term <- replace_terms(res$term, term.replace)
    if (all(is.na(res$nes))) stop("The results have no NES; bars need GSEA results", call. = FALSE)
    if (!is.null(res$comparison)) {
        if (is.null(comparison)) {
            comparison <- levels(res$comparison)[1]
            if (nlevels(res$comparison) > 1) message("Drawing comparison '", comparison, "'; set `comparison` for another")
        }
        if (!(comparison %in% levels(res$comparison))) stop("Comparison '", comparison, "' not found", call. = FALSE)
        res <- res[res$comparison == comparison, , drop = FALSE]
    }
    res <- res[!is.na(res$nes) & !is.na(res$padj), , drop = FALSE]
    res$significant <- res$padj < padj.cutoff
    pick <- function(d, k) {
        if (!show.ns) d <- d[d$significant, , drop = FALSE]
        utils::head(d[order(d$padj, -abs(d$nes)), , drop = FALSE], k)
    }
    df <- rbind(pick(res[res$nes > 0, , drop = FALSE], n.up), pick(res[res$nes < 0, , drop = FALSE], n.down))
    if (nrow(df) == 0) stop("No significant gene sets to draw; see `show.ns`", call. = FALSE)
    df <- df[order(df$nes), , drop = FALSE]
    # Colour value of the adjusted p-values on the scale of `padj.scale`
    df$padj_colour <- padj_colour_value(df$padj, padj.scale)

    sz <- theme_sizes()
    text_size <- sz$small
    # Names sit on the half of the panel opposite their bar: wrap them to fit
    if (is.null(wrap.width)) wrap.width <- max(10, floor((panel.width / 2 - 0.08) / (0.55 * text_size / 72)))
    df$label <- wrap_labels(df$term, wrap.width)
    cell.height <- label_row_height(df$label, cell.height, text_size)
    df$term_id <- make.unique(df$term_id)
    df$y <- factor(df$term_id, levels = df$term_id)
    lim <- max(abs(df$nes)) * 1.05
    # A small gap between 0 and the name, in data units
    gap <- 0.05 * 2 * lim / panel.width
    df$label_x <- ifelse(df$nes > 0, -gap, gap)
    df$hjust <- ifelse(df$nes > 0, 1, 0)

    if (fill == "direction") {
        ends <- bt_palette("blue-magenta")
        dir_cols <- c(up = ends[length(ends)], down = ends[1])
        if (!is.null(cols)) {
            if (is.null(names(cols))) names(cols) <- c("up", "down")[seq_along(cols)]
            dir_cols[names(cols)] <- cols
        }
        df$direction <- factor(ifelse(df$nes > 0, "up", "down"), levels = c("up", "down"))
        bar_fill <- scale_fill_manual(values = dir_cols, guide = "none")
        bar_aes <- aes(fill = .data$direction, alpha = .data$significant)
    } else {
        fill_cols <- resolve_continuous_palette(cols, FALSE, default = "seq-purple")
        # Linear adjusted p-values: the smallest (most significant) darkest
        if (padj.scale == "linear") fill_cols <- rev(fill_cols)
        limits <- c(0, max(df$padj_colour))
        guide <- padj_colour_guide(limits, padj.scale)
        bar_fill <- scale_fill_gradientn(
            colours = fill_cols,
            limits = limits,
            breaks = guide$breaks,
            labels = guide$labels,
            name = bold_title(guide$title),
            guide = guide_colourbar(order = 1, reverse = guide$reverse)
        )
        bar_aes <- aes(fill = .data$padj_colour, alpha = .data$significant)
    }
    # With direction colours, non-significant bars are faded and a legend
    # says so whenever there are any. With adjusted p-value colours the colour
    # already shows significance; fading as well would hide the palest bars
    sig_lab <- significance_labels(padj.cutoff)
    df$significant <- factor(ifelse(df$significant, sig_lab[1], sig_lab[2]), levels = sig_lab)
    fade <- fill == "direction" && any(df$significant == sig_lab[2])
    alpha_scale <- scale_alpha_manual(
        values = stats::setNames(c(1, if (fill == "direction") 0.35 else 1), sig_lab),
        drop = TRUE, name = "Significance",
        guide = if (fade) guide_legend(order = 2, override.aes = list(fill = "grey45")) else "none"
    )

    p <- ggplot(df, aes(x = .data$nes, y = .data$y)) +
        geom_col(bar_aes, width = 0.75) +
        geom_vline(xintercept = 0, colour = "black", linewidth = pt_to_linewidth(sz$line)) +
        geom_text(aes(x = .data$label_x, label = .data$label, hjust = .data$hjust),
                  size = pt_to_size(text_size), lineheight = 0.85, colour = "black") +
        bar_fill +
        alpha_scale +
        scale_x_continuous(limits = c(-lim, lim), expand = expansion(0)) +
        scale_y_discrete(expand = grid_expansion()) +
        labs(x = "Normalized enrichment score", y = NULL) +
        theme_bioToolkit("axis", bg = bg) +
        panel_size_theme(panel.width, grid_panel_size(nrow(df), cell.height)) +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            legend.position = legend.position,
            legend.title = element_text(size = sz$small)
        )
    p
}
