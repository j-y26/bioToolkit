# Shared reader of differential expression results (roadmap item F5)
#
# Design
# ------
# Volcano, MA, DE dot and quadrant plots all start from DE results. Instead of
# each plot guessing column names, `standardize_de()` turns the output of the
# common tools into one table with fixed column names:
#
#   feature, logfc, pvalue, padj, [mean_expr], [comparison], [group], ...
#
# followed by the table's other columns (e.g. Seurat's pct.1 / pct.2).
#
# * Methods are recognised from their column names (presets in
#   R/constants.R); any column can be given explicitly instead, e.g. a shrunken
#   logFC column from DESeq2's lfcShrink().
# * Several comparisons are read from a named list of tables, or from one long
#   table with a comparison column. An optional `groups` mapping adds a second
#   level (e.g. comparison "B cells: KO vs WT" -> group "B cells").
# * p-values of exactly 0 (underflow) are replaced by the smallest non-zero
#   value x 0.1, so that they sit just above the rest of the points on a
#   -log10 axis instead of far off at 1e-300.
# * NA p-values (DESeq2 independent filtering and outlier flagging) are kept
#   as NA and reported; plots treat them as not significant.


#' @title Standardize differential expression results
#'
#' @description Reads DE results from edgeR, DESeq2, limma or Seurat (or any
#'              table, with the columns given explicitly) into one table with
#'              standard column names, so that every DE plot of the package
#'              reads them the same way.
#'
#' @param data A data frame of DE results (or an object coercible to one, such
#'             as a DESeq2 `DESeqResults` or edgeR `TopTags`); or a named list
#'             of them, one per comparison, the names becoming the
#'             `comparison` column
#'
#' @param method The tool that produced the results: "auto" (default) detects
#'               it from the column names and says which was picked; or one of
#'               "edgeR", "DESeq2", "limma" and "seurat"
#'
#' @param feature.col Column of feature (gene) names. If NULL (default), the
#'                    `gene` column of Seurat results when there is one (as
#'                    from \code{FindAllMarkers()}), else the row names
#'
#' @param logfc.col,p.col,padj.col,mean.col Columns of the log2 fold change,
#'        raw p-value, adjusted p-value and mean expression. Each NULL
#'        (default) is taken from `method`. Give e.g. the column of a shrunken
#'        log fold change (DESeq2 \code{lfcShrink()}, apeglm) as `logfc.col`
#'
#' @param comparison.col For a single table holding several comparisons, the
#'                       column naming the comparison of each row (e.g.
#'                       "cluster" for \code{FindAllMarkers()}). Default is
#'                       NULL
#'
#' @param groups Optional second level of comparisons, used for brackets in
#'               multi-comparison plots: a data frame with `comparison` and
#'               `group` columns, or a character vector of groups named by
#'               comparison
#'
#' @details
#' Columns read for each method:
#' \tabular{lllll}{
#'   method \tab logfc \tab pvalue \tab padj \tab mean_expr \cr
#'   edgeR \tab logFC \tab PValue \tab FDR \tab logCPM (log2) \cr
#'   DESeq2 \tab log2FoldChange \tab pvalue \tab padj \tab baseMean (linear) \cr
#'   limma \tab logFC \tab P.Value \tab adj.P.Val \tab AveExpr (log2) \cr
#'   seurat \tab avg_log2FC \tab p_val \tab p_val_adj \tab none
#' }
#'
#' When there is no adjusted p-value column, adjusted p-values are computed
#' with the Benjamini-Hochberg method, per comparison, with a message.
#'
#' p-values of exactly 0, which come from numerical underflow, are replaced by
#' the smallest non-zero p-value of the table times 0.1, with a message, so
#' that they are drawn just beyond the other points rather than at an
#' arbitrary 1e-300. NA p-values (e.g. from DESeq2's independent filtering)
#' are kept and reported; plots treat them as not significant.
#'
#' Seurat versions before 4 reported natural-log fold changes (`avg_logFC`);
#' these are not read automatically, since mixing them with log2 values would
#' be wrong. Convert them (divide by log(2)) and pass the column as
#' `logfc.col`.
#'
#' @return A data frame with the columns `feature`, `logfc`, `pvalue`, `padj`,
#'         `mean_expr` (when available), `comparison` and `group` (when
#'         given), then the other columns of the input. Attributes `method`
#'         (one per comparison) and `mean.scale` ("log2" or "linear", NA when
#'         unknown) record how it was read
#'
#' @examples
#' de <- data.frame(
#'     logFC = c(2.1, -1.4, 0.2), logCPM = c(5, 7, 3),
#'     PValue = c(1e-8, 1e-4, 0.6), FDR = c(3e-8, 2e-4, 0.6),
#'     row.names = c("Cd3e", "Ms4a1", "Actb")
#' )
#' standardize_de(de)
#' standardize_de(list("KO vs WT" = de, "Het vs WT" = de))
#'
#' @export
#'
standardize_de <- function(data,
                           method = c("auto", "edgeR", "DESeq2", "limma", "seurat"),
                           feature.col = NULL,
                           logfc.col = NULL,
                           p.col = NULL,
                           padj.col = NULL,
                           mean.col = NULL,
                           comparison.col = NULL,
                           groups = NULL) {
    method <- match.arg(method)
    cols <- list(feature = feature.col, logfc = logfc.col, pvalue = p.col,
                 padj = padj.col, mean_expr = mean.col)

    if (is.list(data) && !is.data.frame(data) && !inherits(data, "TopTags")) {
        # A list of tables, one per comparison
        if (length(data) == 0) stop("`data` is an empty list", call. = FALSE)
        if (is.null(names(data)) || any(!nzchar(names(data))) || anyDuplicated(names(data))) {
            stop("A list of DE tables must have unique names, which become the comparisons",
                 call. = FALSE)
        }
        if (!is.null(comparison.col)) {
            stop("Give either a named list of tables or `comparison.col`, not both", call. = FALSE)
        }
        parts <- lapply(names(data), function(nm) {
            out <- read_de_table(data[[nm]], method, cols, label = nm)
            out$table$comparison <- nm
            out
        })
        tables <- lapply(parts, `[[`, "table")
        # Columns present in only some of the tables are filled with NA
        all_cols <- unique(unlist(lapply(tables, names)))
        tables <- lapply(tables, function(t) {
            for (col in setdiff(all_cols, names(t))) t[[col]] <- NA
            t[all_cols]
        })
        res <- do.call(rbind, tables)
        res$comparison <- factor(res$comparison, levels = names(data))
        methods <- stats::setNames(vapply(parts, `[[`, "", "method"), names(data))
        # Tables without mean expression (e.g. Seurat) do not decide its scale
        scales_ <- unique(stats::na.omit(vapply(parts, `[[`, "", "mean_scale")))
    } else {
        out <- read_de_table(data, method, cols, label = NULL, comparison.col = comparison.col)
        res <- out$table
        methods <- out$method
        scales_ <- out$mean_scale
    }
    rownames(res) <- NULL

    # Adjusted p-values missing from the input are computed per comparison, so
    # that each comparison is corrected for its own number of tests
    by <- if ("comparison" %in% names(res)) as.character(res$comparison) else rep("", nrow(res))
    computed <- character(0)
    for (b in unique(by)) {
        idx <- which(by == b)
        if (all(is.na(res$padj[idx])) && any(!is.na(res$pvalue[idx]))) {
            res$padj[idx] <- stats::p.adjust(res$pvalue[idx], "BH")
            computed <- c(computed, b)
        }
    }
    if (length(computed)) {
        message("No adjusted p-values", if (any(nzchar(computed))) paste0(" for ", paste(computed, collapse = ", ")),
                ": computed with the Benjamini-Hochberg method")
    }
    res$pvalue <- replace_zero_p(res$pvalue, "p-values")
    res$padj <- replace_zero_p(res$padj, "adjusted p-values")
    n_na <- sum(is.na(res$padj))
    if (n_na > 0) {
        message(n_na, " feature(s) have no adjusted p-value (e.g. DESeq2 independent filtering); ",
                "they are treated as not significant")
    }

    if (!is.null(groups)) res <- add_de_groups(res, groups)

    lead <- intersect(c("feature", "logfc", "pvalue", "padj", "mean_expr", "comparison", "group"), names(res))
    res <- res[c(lead, setdiff(names(res), lead))]
    if (length(scales_) > 1) {
        warning("Mean expression is on different scales across tables (",
                paste(scales_, collapse = ", "), "); it is not comparable between them", call. = FALSE)
    }
    attr(res, "method") <- methods
    attr(res, "mean.scale") <- if (length(scales_) == 1) unname(scales_) else NA_character_
    res
}


#' Read one DE table into standard columns
#'
#' @param data A DE table
#'
#' @param method Method, or "auto"
#'
#' @param cols Columns given explicitly, as a list named by standard column
#'
#' @param label Name of the comparison, for messages, or NULL
#'
#' @param comparison.col Column naming the comparison of each row, or NULL
#'
#' @return A list with the standardized `table`, the `method` used and the
#'         `mean_scale`
#'
#' @keywords internal
read_de_table <- function(data, method, cols, label = NULL, comparison.col = NULL) {
    where <- if (is.null(label)) "" else paste0(" in '", label, "'")
    if (inherits(data, "TopTags")) data <- data$table
    if (!is.data.frame(data)) {
        data <- tryCatch(as.data.frame(data), error = function(e) {
            stop("DE results", where, " must be a data frame", call. = FALSE)
        })
    }
    if (nrow(data) == 0) stop("The DE table", where, " has no rows", call. = FALSE)

    if (method == "auto") method <- detect_de_method(data, cols, where)
    preset <- if (method == "none") {
        c(logfc = NA, pvalue = NA, padj = NA, mean_expr = NA, mean_scale = NA)
    } else {
        DE_COLUMN_PRESETS[[method]]
    }
    pick <- function(std) {
        given <- cols[[std]]
        col <- if (is.null(given)) unname(preset[[std]]) else given
        if (!is.null(given) && !(given %in% names(data))) {
            stop("Column '", given, "' not found", where, call. = FALSE)
        }
        if (is.na(col) || !(col %in% names(data))) NA_character_ else col
    }
    used <- vapply(c("logfc", "pvalue", "padj", "mean_expr"), pick, character(1))
    if (is.na(used[["logfc"]])) {
        stop("No log fold change column", where, "; set `logfc.col`", call. = FALSE)
    }
    if (is.na(used[["pvalue"]]) && is.na(used[["padj"]])) {
        stop("No p-value column", where, "; set `p.col` and / or `padj.col`", call. = FALSE)
    }

    # Feature names: an explicit column, Seurat's gene column (row names of
    # FindAllMarkers() output are made unique, e.g. "Cd3e1"), or the row names
    feature_col <- cols$feature
    if (is.null(feature_col) && method == "seurat" && "gene" %in% names(data)) feature_col <- "gene"
    feature <- if (!is.null(feature_col)) {
        if (!(feature_col %in% names(data))) stop("Column '", feature_col, "' not found", where, call. = FALSE)
        as.character(data[[feature_col]])
    } else {
        # Automatic row names (1, 2, ...) are not feature names
        if (.row_names_info(data) < 0) {
            stop("The DE table", where, " has no row names; set `feature.col`", call. = FALSE)
        }
        rownames(data)
    }

    num <- function(std) {
        col <- used[[std]]
        if (is.na(col)) return(rep(NA_real_, nrow(data)))
        v <- data[[col]]
        if (!is.numeric(v)) stop("Column '", col, "'", where, " is not numeric", call. = FALSE)
        as.numeric(v)
    }
    std <- data.frame(feature = feature, logfc = num("logfc"), pvalue = num("pvalue"),
                      padj = num("padj"), stringsAsFactors = FALSE)
    if (!is.na(used[["mean_expr"]])) std$mean_expr <- num("mean_expr")
    if (is.na(used[["pvalue"]])) message("No raw p-value column", where, "; `pvalue` is NA")

    if (!is.null(comparison.col)) {
        if (!(comparison.col %in% names(data))) {
            stop("Column '", comparison.col, "' not found", where, call. = FALSE)
        }
        cmp <- data[[comparison.col]]
        std$comparison <- if (is.factor(cmp)) droplevels(cmp) else factor(cmp, levels = unique(cmp))
    }
    by <- if (is.null(std$comparison)) rep("", nrow(std)) else std$comparison
    if (anyDuplicated(paste(by, std$feature))) {
        warning("Some features appear more than once", where,
                if (!is.null(comparison.col)) " within a comparison", call. = FALSE)
    }

    # The other columns of the input are kept after the standard ones
    consumed <- c(stats::na.omit(used), feature_col, comparison.col)
    rest <- data[setdiff(names(data), consumed)]
    rest <- rest[setdiff(names(rest), names(std))]
    table <- if (ncol(rest)) cbind(std, rest, stringsAsFactors = FALSE) else std

    mean_scale <- if (!is.null(cols$mean_expr) || is.na(used[["mean_expr"]]) || method == "none") {
        NA_character_
    } else {
        unname(preset[["mean_scale"]])
    }
    list(table = table, method = method, mean_scale = mean_scale)
}


#' Detect the tool that produced a DE table from its column names
#'
#' @param data A DE table
#'
#' @param cols Columns given explicitly
#'
#' @param where Location for messages
#'
#' @return A method name, or "none" when the columns were all given explicitly
#'
#' @keywords internal
detect_de_method <- function(data, cols, where) {
    found <- names(Filter(function(p) {
        all(c(p[["logfc"]], p[["pvalue"]], p[["padj"]]) %in% names(data))
    }, DE_COLUMN_PRESETS))
    if (length(found) >= 1) {
        message("Detected ", found[1], " results", where, " (",
                paste(DE_COLUMN_PRESETS[[found[1]]][c("logfc", "pvalue", "padj")], collapse = ", "), ")")
        return(found[1])
    }
    # Columns given explicitly need no method
    if (!is.null(cols$logfc) && (!is.null(cols$pvalue) || !is.null(cols$padj))) return("none")
    if ("avg_logFC" %in% names(data)) {
        stop("Column 'avg_logFC'", where, " holds natural-log fold changes from Seurat < 4; ",
             "convert it to log2 (divide by log(2)) and pass it as `logfc.col`", call. = FALSE)
    }
    stop("Could not detect the DE method", where, " from its columns (",
         paste(utils::head(names(data), 8), collapse = ", "),
         "); set `method` or the column arguments", call. = FALSE)
}


#' Replace p-values of exactly 0 by a tenth of the smallest non-zero one
#'
#' @param p p-values
#'
#' @param what Name of the values, for the message
#'
#' @return The p-values
#'
#' @keywords internal
replace_zero_p <- function(p, what) {
    zero <- !is.na(p) & p == 0
    if (!any(zero)) return(p)
    positive <- p[!is.na(p) & p > 0]
    floor_p <- if (length(positive)) min(positive) * 0.1 else .Machine$double.xmin
    message(sum(zero), " ", what, " of 0 replaced by ", signif(floor_p, 3),
            " (a tenth of the smallest non-zero value)")
    p[zero] <- floor_p
    p
}


#' Add a group to each comparison
#'
#' @param res Standardized DE results with a comparison column
#'
#' @param groups A data frame with comparison and group columns, or a
#'               character vector of groups named by comparison
#'
#' @return `res` with a `group` factor, its levels in the order the groups
#'         first appear along the comparisons
#'
#' @keywords internal
add_de_groups <- function(res, groups) {
    if (is.null(res$comparison)) {
        stop("`groups` needs several comparisons (a list of tables or `comparison.col`)", call. = FALSE)
    }
    if (is.data.frame(groups)) {
        if (!all(c("comparison", "group") %in% names(groups))) {
            stop("`groups` must have `comparison` and `group` columns", call. = FALSE)
        }
        groups <- stats::setNames(as.character(groups$group), as.character(groups$comparison))
    }
    if (is.null(names(groups))) stop("`groups` must be named by comparison", call. = FALSE)
    cmp <- levels(res$comparison)
    missing <- setdiff(cmp, names(groups))
    if (length(missing)) {
        stop("No group given for comparison(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
    group_of <- groups[cmp]
    res$group <- factor(unname(groups[as.character(res$comparison)]), levels = unique(group_of))
    res
}
