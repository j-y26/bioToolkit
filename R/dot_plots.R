# Marker dot plot for Seurat objects (roadmap item G1)
#
# Design
# ------
# `DotPlot_bioToolkit()` replaces Seurat::DotPlot() and covers the features of
# SCpubr::do_DotPlot(): mean and scaled (z-scored) mean expression, gene
# groups from a named list, `split.by`, clustering of genes and groups,
# `flip`. It is drawn by the dot-grid engine (G0, R/dot_grid.R).
#
# * Assay: read explicitly and reported (R/seurat_data.R). Seurat and SCpubr
#   both read the default assay; with "RNA" and "SCT" in one object the same
#   call then gives different numbers depending on `DefaultAssay()`.
# * Mean expression: two definitions are in use (`average`):
#     "linear" (default, as Seurat::DotPlot): log1p(mean(expm1(x))), the log of
#              the mean normalized count. It is the expression of the group as
#              a whole (a pseudobulk), on the same scale as the fold changes
#              of Seurat::FindMarkers() (avg_log2FC is computed this way).
#     "log"    (as SCpubr::do_DotPlot): mean(x), the mean of the
#              log-normalized values. Less driven by a few cells with very
#              high counts; closer to the "typical" cell.
#   Both only apply to the "data" layer (log1p values); other layers are
#   averaged as they are.
# * Scaling (`scale = TRUE`, default as in Seurat): each gene's means are
#   z-scored across the plotted groups (and splits), which compares groups
#   for each gene but hides how highly genes are expressed overall. The colour
#   scale is diverging, symmetric around 0, and reaches its ends at
#   +-min(max |z|, 2.5) (Seurat clips at +-2.5); source data keep the
#   unclipped values.
# * % expressing: the share of cells with a value above `expr.threshold`
#   (default 0), as Seurat. Dot areas are proportional to it, from 0 to
#   `pct.max` (default 100%), so that a dot of a given size means the same
#   percentage in every figure (Seurat and SCpubr rescale sizes to each plot's
#   own range).
# * Colours: scaled means on the diverging "blue-magenta"; unscaled means on
#   "viridis_r" (yellow = low, dark purple = high; SCpubr's viridis option
#   runs the same way).
# * Groups: a marker next to each group label, in the default group colours
#   (as DimPlot_bioToolkit()) or `group.cols`, so that the dot plot can be
#   matched with an embedding by colour: a coloured dot (default) or a
#   colour strip (`group.marker`).


#' @title Marker dot plot
#'
#' @description Shows, for each gene and group of cells, the percentage of
#'              cells expressing the gene (dot size) and its mean expression
#'              (colour), scaled across groups by default. A drop-in
#'              alternative to \code{Seurat::DotPlot()} and
#'              \code{SCpubr::do_DotPlot()}, with gene groups, clustering, a
#'              coloured marker for each group, gene names in italics and a fixed
#'              size: every row and column is `cell.height` /
#'              `cell.width` inches, so dots and text are the same size in
#'              every figure. Save it with \code{\link{save_plot}}; the
#'              numbers behind it are returned by
#'              \code{\link{export_source_data}}.
#'
#' @param object A Seurat object
#'
#' @param features Genes to plot: a character vector, or a named list of
#'                 character vectors (e.g. markers per cell type, see
#'                 \code{\link{top_markers}}), drawn as labelled groups of
#'                 columns. Genes are matched exactly (rat and mouse symbols
#'                 are e.g. "Cd3e", human "CD3E"); genes not in the assay are
#'                 left out with a warning. A gene in two groups is drawn in
#'                 both
#'
#' @param group.by Metadata column to group cells by. If NULL (default), the
#'                 active identities
#'
#' @param split.by Metadata column to split each group by (e.g. condition),
#'                 giving one row per group and split, with the rows of a group
#'                 together. Default is NULL
#'
#' @param cells Cells to use, as names or indices. If NULL (default), all
#'
#' @param assay Assay to read expression from. If NULL (default), the default
#'              assay of the object, as in Seurat and SCpubr. The assay and
#'              layer used are always reported in a message and written to the
#'              source data. Set it explicitly when the object holds several
#'              assays (e.g. "RNA" and "SCT"), which give different values for
#'              the same gene
#'
#' @param slot Layer to read, default is "data" (log-normalized values). For
#'             "data", see `average`; other layers are averaged as they are
#'
#' @param average How the mean expression of a group is computed from the
#'                log-normalized "data" layer:
#'                \describe{
#'                  \item{"linear"}{(default, as \code{Seurat::DotPlot()})
#'                        \code{log1p(mean(expm1(x)))}: the log of the mean
#'                        normalized expression, i.e. a pseudobulk of the
#'                        group, on the scale of Seurat's fold changes}
#'                  \item{"log"}{(as \code{SCpubr::do_DotPlot()})
#'                        \code{mean(x)}: the mean of the log-normalized
#'                        values, less influenced by a few cells with very high
#'                        expression}
#'                }
#'
#' @param scale Whether to z-score each gene's mean expression across the
#'              plotted groups, default is TRUE (as Seurat; SCpubr's
#'              `zscore.data`). Scaled colours compare groups for each gene,
#'              not genes with each other. Needs at least two groups
#'
#' @param expr.threshold A cell expresses a gene when its value is above this,
#'                       default is 0
#'
#' @param pct.min Dots of groups with fewer than this percentage of cells
#'                expressing the gene are not drawn (SCpubr's `dot.min`, which
#'                defaults to 5), default is 0. They still count for scaling
#'                and are kept in the source data
#'
#' @param pct.max Percentage drawn as the largest dot, default is 100. Lower
#'                it for genes detected in few cells (e.g. transcription
#'                factors); higher percentages are drawn at the largest size
#'                and the legend says so
#'
#' @param min.cutoff,max.cutoff Ends of the colour scale: a number, a quantile
#'                              of the plotted values such as "q5" / "q95",
#'                              or NA (default). With `scale = TRUE`, NA gives
#'                              a scale symmetric around 0 reaching
#'                              +-min(largest |z-score|, 2.5); otherwise 0
#'                              and the largest mean. Values beyond the ends
#'                              get the colour of the end
#'
#' @param cluster.features Whether to order genes by hierarchical clustering
#'                         of their colour values (within each group of
#'                         `features`), default is FALSE
#'
#' @param cluster.groups Whether to order the groups by hierarchical
#'                       clustering of their colour values, default is FALSE.
#'                       With `split.by`, the rows of a group stay together
#'
#' @param cols Colours of the expression scale. If NULL (default),
#'             "blue-magenta" (diverging) with `scale = TRUE`, else
#'             "viridis_r" (sequential, yellow to dark purple). Otherwise a
#'             palette name (see \code{\link{show_palettes}}) or colours
#'             from low to high, e.g. `c("grey90", "darkred")` or
#'             `c(low, mid, high)`
#'
#' @param group.marker How each group is marked with its colour next to its
#'                     label: "dot" (default; a coloured dot), "strip" (a
#'                     colour strip) or "none"
#'
#' @param group.cols Colours of the group markers: NULL (default) for
#'                   the colours \code{\link{DimPlot_bioToolkit}} gives them
#'                   ("scanpy-20" up to 20 groups), a qualitative palette name,
#'                   or colours named by group (or in level order)
#'
#' @param flip Whether to draw genes as rows and groups as columns, default is
#'             FALSE (genes as columns, as Seurat)
#'
#' @param italic.genes Whether gene names are set in italics, default is TRUE.
#'                     Set FALSE for features that are not genes (e.g.
#'                     proteins of an ADT assay)
#'
#' @param fill.title Title of the colour legend. If NULL (default), "Mean
#'                   expression", with "(z-score)" when scaled
#'
#' @param size.title Title of the size legend, default is "\% expressing"
#'
#' @param cell.width,cell.height Width of a column and height of a row in
#'                               inches, default is 0.18 (13 pt). Dots scale
#'                               with them
#'
#' @param legend.position Position of the legends: "right" (default),
#'                        "left", "top", "bottom" or "none"
#'
#' @param bg Background colour, default is "transparent"
#'
#' @details
#' Cells without a value of `group.by` (or `split.by`) are left out, with a
#' message. Groups (and group-split combinations) with no cells are not
#' drawn; group colours are still taken from all levels of `group.by` in the
#' whole object, so that they match an embedding of all cells.
#'
#' The source data (\code{\link{export_source_data}}) hold one row per group
#' (and split) and gene: the number of cells, `pct_expressing`,
#' `mean_expression`, `mean_expression_scaled` (when scaled), whether the dot
#' is `drawn`, and the `assay`, `layer` and `average` used.
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' # Genes as a vector, from the SCT assay
#' DotPlot_bioToolkit(obj, features = c("Cd3e", "Ms4a1", "Cd68"),
#'                    group.by = "cell_type", assay = "SCT")
#'
#' # Marker sets per cell type, drawn as labelled groups of columns
#' markers <- list(
#'     "T cells" = c("Cd3e", "Cd3d"),
#'     "B cells" = c("Ms4a1", "Cd79a"),
#'     "Macrophages" = c("Cd68", "Lyz2")
#' )
#' p <- DotPlot_bioToolkit(obj, features = markers, group.by = "cell_type")
#' save_plot(p, "dotplot.pdf")
#'
#' # As SCpubr::do_DotPlot(): mean of log values, unscaled, dots from 5\%
#' DotPlot_bioToolkit(obj, features = markers, group.by = "cell_type",
#'                    average = "log", scale = FALSE, pct.min = 5)
#'
#' # Top markers from FindAllMarkers()
#' all_markers <- Seurat::FindAllMarkers(obj, only.pos = TRUE)
#' DotPlot_bioToolkit(obj, features = top_markers(all_markers, n = 3))
#' }
#'
#' @seealso \code{\link{top_markers}}, \code{\link{FeaturePlot_bioToolkit}},
#'          \code{\link{save_plot}}, \code{\link{export_source_data}}
#'
#' @import ggplot2
#' @export
#'
DotPlot_bioToolkit <- function(
    object,
    features,
    group.by = NULL,
    split.by = NULL,
    cells = NULL,
    assay = NULL,
    slot = "data",
    average = c("linear", "log"),
    scale = TRUE,
    expr.threshold = 0,
    pct.min = 0,
    pct.max = 100,
    min.cutoff = NA,
    max.cutoff = NA,
    cluster.features = FALSE,
    cluster.groups = FALSE,
    cols = NULL,
    group.marker = c("dot", "strip", "none"),
    group.cols = NULL,
    flip = FALSE,
    italic.genes = TRUE,
    fill.title = NULL,
    size.title = "% expressing",
    cell.width = 0.18,
    cell.height = 0.18,
    legend.position = c("right", "left", "top", "bottom", "none"),
    bg = "transparent"
) {
    check_installed("Seurat", reason = "to fetch expression")
    average <- match.arg(average)
    group.marker <- match.arg(group.marker)
    legend.position <- match.arg(legend.position)
    if (!(is.numeric(pct.max) && length(pct.max) == 1 && pct.max > 0 && pct.max <= 100)) {
        stop("`pct.max` must be a single number in (0, 100]", call. = FALSE)
    }
    if (!(is.numeric(pct.min) && length(pct.min) == 1 && pct.min >= 0 && pct.min < 100)) {
        stop("`pct.min` must be a single number in [0, 100)", call. = FALSE)
    }

    # Genes, with their groups when `features` is a named list
    ft <- dot_feature_table(features)
    assay_info <- resolve_assay(object, assay, slot)
    genes <- genes_in_assay(object, unique(ft$feature), assay_info)
    if (length(genes) == 0) stop("None of `features` are in assay '", assay_info$assay, "'", call. = FALSE)
    ft <- ft[ft$feature %in% genes, , drop = FALSE]

    # Cells and their groups. Group levels come from the whole object, so
    # that strip colours match an embedding of all cells
    cells <- resolve_cells(cells, colnames(object))
    if (length(cells) == 0) stop("No cells to plot", call. = FALSE)
    group_all <- fetch_group(object, group.by, colnames(object))
    group <- group_all[match(cells, colnames(object))]
    split <- NULL
    if (!is.null(split.by)) {
        if (length(split.by) != 1) stop("`split.by` must be a single metadata column", call. = FALSE)
        split <- fetch_meta(object, split.by, cells)
        if (!is.factor(split)) split <- factor(split)
    }
    keep <- !is.na(group) & (if (is.null(split)) TRUE else !is.na(split))
    if (!all(keep)) {
        message(sum(!keep), " cell(s) without a value of `",
                paste(c(if (is.null(group.by)) "identity" else group.by, split.by), collapse = "` or `"),
                "` are left out")
        cells <- cells[keep]
        group <- group[keep]
        if (!is.null(split)) split <- split[keep]
    }

    # One unit (row of the plot by default) per group, or per group and split
    units <- dot_units(group, split)
    expr <- as.matrix(fetch_features(object, genes, cells, slot, assay_info))
    if (any(expr < 0, na.rm = TRUE)) {
        warning("Layer '", slot, "' of assay '", assay_info$assay, "' has negative values: it is not ",
                "normalized expression (e.g. scaled or integrated data), so the % of cells expressing and ",
                "the means are not meaningful. Use the 'data' layer of a normalized assay (e.g. RNA or SCT)",
                call. = FALSE)
    }
    st <- group_expression_stats(expr, units$unit, average, expr.threshold, log.data = slot == "data")

    if (scale && nrow(st$mean) < 2) {
        message("Only one group: means are not scaled")
        scale <- FALSE
    }
    fill_mat <- if (scale) scale_columns(st$mean) else st$mean

    # Long table: one row per unit and gene (a gene in two groups twice)
    df <- dot_long_table(ft, units$table, st, fill_mat, scale)
    df$drawn <- df$pct_expressing >= pct.min

    # Order of genes (within groups) and of groups
    feat_lv <- unique(ft$id)
    if (cluster.features) {
        m <- t(fill_mat[, ft$feature, drop = FALSE])
        rownames(m) <- ft$id
        feat_lv <- cluster_order(m, ft$feature_group)
    }
    group_lv <- levels(droplevels(units$table$group))
    if (cluster.groups) {
        # One profile per group: its values for every gene (and split)
        key <- if (is.null(df$split)) df$id else paste(df$id, df$split)
        m <- tapply(df$fill, list(as.character(df$group), key), mean)
        group_lv <- cluster_order(m[group_lv, , drop = FALSE])
    }
    unit_tab <- units$table
    unit_tab$group <- factor(as.character(unit_tab$group), levels = group_lv)
    ord <- if (is.null(unit_tab$split)) order(unit_tab$group) else order(unit_tab$group, unit_tab$split)
    unit_tab <- unit_tab[ord, , drop = FALSE]
    unit_lv <- as.character(unit_tab$unit)

    # Colour limits
    limits <- dot_fill_limits(df$fill[df$drawn], scale, min.cutoff, max.cutoff)
    # Size scale: areas from 0 to pct.max, legend labels with ">=" at the top
    # when larger percentages are squished onto it
    size_breaks <- scales::breaks_extended(5)(c(0, pct.max))
    size_breaks <- size_breaks[size_breaks > 0 & size_breaks <= pct.max]
    size_labels <- format(size_breaks, trim = TRUE)
    if (any(df$pct_expressing[df$drawn] > pct.max)) {
        if (utils::tail(size_breaks, 1) != pct.max) {
            size_breaks <- c(size_breaks, pct.max)
            size_labels <- c(size_labels, format(pct.max, trim = TRUE))
        }
        size_labels[length(size_labels)] <- paste0("\u2265", size_labels[length(size_labels)])
    }

    # Group colours for the markers
    has_marker <- group.marker != "none"
    strip_cols <- if (has_marker) group_colours(group.cols, droplevels(group_all), arg = "group.cols")$cols

    # Engine table: genes as columns (default) or rows (flip)
    plot_df <- df
    feat <- factor(plot_df$id, levels = feat_lv)
    unit <- factor(as.character(plot_df$unit), levels = unit_lv)
    feat_group <- if (!is.null(ft$feature_group)) factor(plot_df$feature_group, levels = levels(ft$feature_group))
    unit_group <- if (!is.null(split)) factor(as.character(plot_df$group), levels = group_lv)
    feat_labels <- stats::setNames(ft$feature, ft$id)
    unit_labels <- stats::setNames(
        as.character(if (is.null(split)) unit_tab$group else unit_tab$split), unit_lv
    )
    # Dots below `pct.min` are not drawn (NA size), but their rows and
    # columns stay
    eng <- data.frame(size = ifelse(plot_df$drawn, plot_df$pct_expressing, NA), fill = plot_df$fill)
    if (!flip) {
        eng$x <- feat
        eng$y <- unit
        eng$x_group <- feat_group
        eng$y_group <- unit_group
    } else {
        eng$x <- unit
        eng$y <- feat
        eng$x_group <- unit_group
        eng$y_group <- feat_group
    }
    if (has_marker) eng$strip <- unname(strip_cols[as.character(plot_df$group)])
    if (is.null(fill.title)) fill.title <- if (scale) "Mean expression\n(z-score)" else "Mean\nexpression"
    italic <- stats::setNames(rep(italic.genes, length(feat_labels)), names(feat_labels))

    p <- dot_grid_plot(
        eng,
        x.labels = if (!flip) feat_labels else unit_labels,
        y.labels = if (!flip) unit_labels else feat_labels,
        x.italic = if (!flip) italic else FALSE,
        y.italic = if (!flip) FALSE else italic,
        size.limits = c(0, pct.max),
        size.breaks = size_breaks,
        size.labels = size_labels,
        size.title = size.title,
        fill.colours = resolve_continuous_palette(cols, diverging = scale, default = "viridis_r"),
        fill.limits = limits,
        fill.center = if (scale) 0 else NULL,
        fill.title = fill.title,
        strip.axis = if (has_marker) (if (!flip) "y" else "x"),
        strip.style = if (group.marker == "dot") "dot" else "strip",
        cell.width = cell.width,
        cell.height = cell.height,
        legend.position = legend.position,
        bg = bg
    )

    # Source data: every unit and gene, drawn or not, with how it was computed
    src <- df[, intersect(c("group", "split", "feature_group", "feature", "n_cells", "pct_expressing",
                            "mean_expression", "mean_expression_scaled", "drawn"), names(df))]
    src$assay <- assay_info$assay
    src$layer <- slot
    src$average <- if (slot == "data") average else "arithmetic"
    attr(p, "source_data") <- src
    p
}


#' Genes of a dot plot, with their groups
#'
#' @param features A character vector, or a named list of character vectors
#'
#' @return A data frame with `feature`, `id` (unique per group and gene) and,
#'         for a list, `feature_group` (a factor in list order)
#'
#' @keywords internal
dot_feature_table <- function(features) {
    if (is.list(features)) {
        if (is.null(names(features)) || any(!nzchar(names(features))) || anyDuplicated(names(features))) {
            stop("A list of `features` must have unique names, which label the groups", call. = FALSE)
        }
        features <- lapply(features, function(f) unique(as.character(unlist(f))))
        features <- features[lengths(features) > 0]
        if (length(features) == 0) stop("No `features` given", call. = FALSE)
        ft <- data.frame(
            feature = unlist(features, use.names = FALSE),
            feature_group = factor(rep(names(features), lengths(features)), levels = names(features)),
            stringsAsFactors = FALSE
        )
        dup <- unique(ft$feature[duplicated(ft$feature)])
        if (length(dup)) message("Gene(s) in more than one group, drawn in each: ", paste(dup, collapse = ", "))
        ft$id <- paste(ft$feature_group, ft$feature, sep = "\r")
        return(ft)
    }
    features <- unique(as.character(features))
    if (length(features) == 0) stop("No `features` given", call. = FALSE)
    data.frame(feature = features, id = features, stringsAsFactors = FALSE)
}


#' Plotting units of a dot plot: groups, or group and split combinations
#'
#' @param group A factor along the cells
#'
#' @param split A factor along the cells, or NULL
#'
#' @return A list with `unit`, a factor along the cells, and `table`, one row
#'         per unit with cells (`unit`, `group` and, when split, `split`), in
#'         group then split level order
#'
#' @keywords internal
dot_units <- function(group,
                      split) {
    group <- droplevels(group)
    if (is.null(split)) {
        return(list(unit = group, table = data.frame(unit = factor(levels(group), levels = levels(group)),
                                                     group = factor(levels(group), levels = levels(group)))))
    }
    split <- droplevels(split)
    combos <- expand.grid(split = levels(split), group = levels(group), stringsAsFactors = FALSE)
    combos$unit <- paste(combos$group, combos$split, sep = "\r")
    unit <- factor(paste(group, split, sep = "\r"), levels = combos$unit)
    present <- combos$unit %in% unit
    combos <- combos[present, , drop = FALSE]
    unit <- droplevels(unit)
    list(
        unit = unit,
        table = data.frame(
            unit = factor(combos$unit, levels = combos$unit),
            group = factor(combos$group, levels = levels(group)),
            split = factor(combos$split, levels = levels(split))
        )
    )
}


#' Mean expression and percentage of expressing cells per group
#'
#' @param expr A matrix of cells x genes
#'
#' @param unit A factor along the cells
#'
#' @param average "linear" or "log"
#'
#' @param threshold A cell expresses a gene above this value
#'
#' @param log.data Whether `expr` holds log1p values (the "data" layer), for
#'                 which `average` applies
#'
#' @return A list of three matrices, units x genes: `mean`, `pct` and `n`
#'         (number of cells)
#'
#' @keywords internal
group_expression_stats <- function(expr,
                                   unit,
                                   average,
                                   threshold,
                                   log.data) {
    n <- as.vector(table(unit))
    linear <- log.data && average == "linear"
    # rowsum() sums the rows of each unit; dividing by the number of cells
    # gives the means
    sums <- rowsum(if (linear) expm1(expr) else expr, unit, reorder = TRUE)
    mean <- sums / n
    if (linear) mean <- log1p(mean)
    pct <- rowsum((expr > threshold) * 1, unit, reorder = TRUE) / n * 100
    list(mean = mean, pct = pct, n = stats::setNames(n, levels(unit)))
}


#' z-score the columns of a matrix
#'
#' @param m A matrix (units x genes)
#'
#' @return The matrix with each column centred and divided by its standard
#'         deviation; constant columns become 0
#'
#' @keywords internal
scale_columns <- function(m) {
    z <- scale(m)
    z[!is.finite(z)] <- 0
    attr(z, "scaled:center") <- NULL
    attr(z, "scaled:scale") <- NULL
    z
}


#' Long table of a dot plot
#'
#' @param ft Gene table from \code{dot_feature_table()}
#'
#' @param units Unit table from \code{dot_units()}
#'
#' @param st Statistics from \code{group_expression_stats()}
#'
#' @param fill_mat Values the dots are coloured by, units x genes
#'
#' @param scale Whether `fill_mat` holds z-scores
#'
#' @return A data frame, one row per unit and row of `ft`
#'
#' @keywords internal
dot_long_table <- function(ft,
                           units,
                           st,
                           fill_mat,
                           scale) {
    u <- as.character(units$unit)
    idx <- expand.grid(i = seq_len(nrow(units)), j = seq_len(nrow(ft)))
    df <- data.frame(unit = u[idx$i], group = units$group[idx$i], stringsAsFactors = FALSE)
    if (!is.null(units$split)) df$split <- units$split[idx$i]
    if (!is.null(ft$feature_group)) df$feature_group <- ft$feature_group[idx$j]
    df$feature <- ft$feature[idx$j]
    df$id <- ft$id[idx$j]
    pick <- function(m) m[cbind(match(df$unit, rownames(m)), match(df$feature, colnames(m)))]
    df$n_cells <- unname(st$n[df$unit])
    df$pct_expressing <- pick(st$pct)
    df$mean_expression <- pick(st$mean)
    if (scale) df$mean_expression_scaled <- pick(fill_mat)
    df$fill <- pick(fill_mat)
    df
}


#' Colour limits of a dot plot
#'
#' @param v Values coloured
#'
#' @param scale Whether they are z-scores
#'
#' @param min.cutoff,max.cutoff Cutoffs (number, "q##" or NA)
#'
#' @return Two numbers
#'
#' @keywords internal
dot_fill_limits <- function(v,
                            scale,
                            min.cutoff,
                            max.cutoff) {
    v <- v[is.finite(v)]
    if (length(v) == 0) return(c(0, 1))
    if (scale) {
        # Symmetric around 0, at most +-2.5 (Seurat's col.min / col.max)
        m <- min(max(abs(v)), 2.5)
        if (m == 0) m <- 1
        default_lo <- function(x) -m
        default_hi <- function(x) m
    } else {
        default_lo <- function(x) min(0, min(x))
        default_hi <- function(x) max(x)
    }
    lim <- c(cutoff_value(min.cutoff, v, default_lo), cutoff_value(max.cutoff, v, default_hi))
    if (lim[2] <= lim[1]) lim[2] <- lim[1] + 1
    lim
}


#' @title Top marker genes per group
#'
#' @description Picks the top marker genes of each group from
#'              \code{Seurat::FindAllMarkers()} output (or any DE table with a
#'              group column), as a named list for the `features` of
#'              \code{\link{DotPlot_bioToolkit}}.
#'
#' @param markers A data frame of markers, e.g. from
#'                \code{Seurat::FindAllMarkers()}; read with
#'                \code{\link{standardize_de}}
#'
#' @param n Number of genes per group, default is 5
#'
#' @param group.col Column naming the group of each row, default is "cluster"
#'                  (as FindAllMarkers)
#'
#' @param rank.by How genes are ranked within a group: "logfc" (default;
#'                largest log2 fold change, ties broken by the adjusted
#'                p-value) or "padj" (smallest adjusted p-value, ties broken by
#'                the largest log2 fold change). Ranking by p-value favours
#'                broadly expressed genes in large groups; ranking by fold
#'                change favours specific genes
#'
#' @param padj.cutoff Only genes with an adjusted p-value below this, and a
#'                    positive log2 fold change, are markers. Default is 0.05
#'
#' @param unique Whether each gene is a marker of one group only, default is
#'               TRUE: a gene that is a marker of several groups is kept for
#'               the group where its log2 fold change is largest
#'
#' @param ... Further arguments to \code{\link{standardize_de}}, e.g.
#'            `feature.col`
#'
#' @return A named list of character vectors, one per group with markers, in
#'         the order of the groups
#'
#' @examples
#' markers <- data.frame(
#'     gene = c("Cd3e", "Cd3d", "Ms4a1", "Cd79a", "Cd3e"),
#'     cluster = factor(c("T", "T", "B", "B", "B")),
#'     avg_log2FC = c(3, 2.5, 4, 3.5, 0.5),
#'     p_val = 1e-10, p_val_adj = 1e-8
#' )
#' top_markers(markers, n = 2)
#'
#' @seealso \code{\link{DotPlot_bioToolkit}}
#'
#' @export
#'
top_markers <- function(markers,
                        n = 5,
                        group.col = "cluster",
                        rank.by = c("logfc", "padj"),
                        padj.cutoff = 0.05,
                        unique = TRUE,
                        ...) {
    rank.by <- match.arg(rank.by)
    if (!(is.numeric(n) && length(n) == 1 && n >= 1)) stop("`n` must be a single number >= 1", call. = FALSE)
    de <- suppressMessages(standardize_de(markers, comparison.col = group.col, ...))
    de <- de[!is.na(de$padj) & de$padj < padj.cutoff & !is.na(de$logfc) & de$logfc > 0, , drop = FALSE]
    if (nrow(de) == 0) stop("No positive markers with an adjusted p-value below ", padj.cutoff, call. = FALSE)
    if (unique) {
        # Each gene stays with the group where it changes most
        de <- de[order(-de$logfc, de$padj), , drop = FALSE]
        de <- de[!duplicated(de$feature), , drop = FALSE]
    }
    ord <- if (rank.by == "logfc") order(-de$logfc, de$padj) else order(de$padj, -de$logfc)
    de <- de[ord, , drop = FALSE]
    out <- lapply(split(de$feature, de$comparison), utils::head, n)
    empty <- lengths(out) == 0
    if (any(empty)) message("No markers for: ", paste(names(out)[empty], collapse = ", "))
    out[!empty]
}
