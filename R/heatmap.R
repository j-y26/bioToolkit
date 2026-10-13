# Marker heatmap for Seurat objects (roadmap item S2)
#
# Design
# ------
# `DoHeatmap_bioToolkit()` replaces Seurat::DoHeatmap(). Decisions (user,
# 2026-10-13):
#
# * Columns are cells (`columns = "cells"`) or group means
#   (`columns = "groups"`), arranged in two levels, as every plot of the
#   package: `split.by` (higher level, e.g. cell type) makes blocks of
#   columns, named on top and set apart by gaps, like facets; `group.by`
#   (lower level, e.g. condition) gives the columns compared within a block,
#   marked by a coloured dot (default), a colour strip or nothing, with a
#   legend, since such columns are too narrow to name. Without `split.by`,
#   the groups are the blocks (named on top) and the markers show the
#   DimPlot_bioToolkit() colours.
# * No clustering and no dendrograms (use ComplexHeatmap for that). Genes are
#   in the order of `features`, in blocks for a named list (e.g. from
#   top_markers()); cells are in random order within their group, or ordered
#   by the first of `annotations`.
# * Values. Seurat::DoHeatmap() reads "scale.data", which holds only the
#   variable genes (others are dropped), is computed on all cells of the
#   object (stale after subsetting) and, for SCTransform, holds scaled
#   Pearson residuals rather than expression. Here the default is the "data"
#   layer (every gene), z-scored per feature (`scale = TRUE`):
#     cells   across all plotted cells, before downsampling, so that colours
#             do not depend on which cells are drawn;
#     groups  across the plotted columns (group means), as DotPlot_bioToolkit().
#   `slot = "scale.data"` is drawn as stored; `scale = FALSE` draws values
#   as they are (module scores, chromVAR deviations). Colours are clipped at
#   +-2.5 (Seurat's disp.min / disp.max) through `min.cutoff` /
#   `max.cutoff`, as in DotPlot_bioToolkit().
# * Colours (user choice, 2026-10-13): for scaled or signed values, centred
#   on 0, "teal-yellow" (teal - black - yellow) for columns of cells, where
#   most tiles are near 0 and a dark centre lets the few high cells stand
#   out; "blue-magenta" (light grey centre, as DotPlot_bioToolkit()) for
#   group means, a few large tiles read best on a light background. "magma"
#   for unscaled non-negative values (dark = low).
# * Downsampling (`max.cells` per group and block, seeded) keeps small groups
#   visible: Seurat::DoHeatmap() gives each cell the same width, so large
#   groups take most of the figure.
# * Gene names are drawn when a row is tall enough for them (>= the 9 pt
#   text, 0.125 in), else hidden with a message; `label.features` names a
#   few genes with leader lines instead.
#
# Geometry (F2). The body is one ggplot of geom_raster() tiles (one image per
# panel in a PDF, not one rectangle per tile), one facet per block of columns
# and per block of genes. Columns and rows sit at integer positions, so that
# the raster grid is regular; groups within a block are set apart by empty
# columns. Header rows (group markers, annotation strips) and the gene
# labels with leader lines are separate small ggplots with the same facets
# and widths, joined with patchwork, so that each has its own colour scale
# and legend.


#' @title Heatmap of genes across cells or groups
#'
#' @description Shows the expression of genes (e.g. top markers of each cell
#'              type) across cells, or across group means, z-scored per gene by
#'              default. Columns are arranged in blocks (`split.by`, e.g. cell
#'              type) and, within them, groups compared with each other
#'              (`group.by`, e.g. condition) marked by a coloured dot or strip.
#'              A drop-in alternative to \code{Seurat::DoHeatmap()} that reads
#'              any gene (not only those in "scale.data"), keeps small groups
#'              visible by downsampling, and has a fixed size. Save it with
#'              \code{\link{save_plot}}; the values behind it are returned by
#'              \code{\link{export_source_data}}.
#'
#' @param object A Seurat object
#'
#' @param features Genes (or other numeric features, e.g. module scores) as
#'                 rows: a character vector, or a named list of character
#'                 vectors (e.g. from \code{\link{top_markers}}), drawn as
#'                 blocks of rows named by the list. Genes not in the assay
#'                 are left out with a warning
#'
#' @param group.by Metadata column of the groups compared (e.g. condition),
#'                 marked by a coloured dot or strip above the columns. If
#'                 NULL (default), the active identities
#'
#' @param split.by Metadata column of a higher-level grouping (e.g. cell
#'                 type), drawn like facets: one named block of columns per
#'                 level. If NULL (default), the groups are the blocks
#'
#' @param cells Cells to use, as names or indices. If NULL (default), all
#'
#' @param assay Assay to read genes from. If NULL (default), the default
#'              assay. The assay and layer used are reported in a message and
#'              written to the source data
#'
#' @param slot Layer to read, default is "data" (log-normalized values, all
#'             genes). "scale.data" is drawn as stored (`scale` is then not
#'             applied); it usually holds only the variable genes
#'
#' @param columns What a column is: "cells" (default; one column per cell,
#'                at most `max.cells` per group in each block) or "groups"
#'                (one column per group in each block: the group's mean)
#'
#' @param average How group means are computed from the log-normalized
#'                "data" layer (`columns = "groups"`): "linear" (default,
#'                \code{log1p(mean(expm1(x)))}, as Seurat) or "log"
#'                (\code{mean(x)}), as in \code{\link{DotPlot_bioToolkit}}
#'
#' @param max.cells Largest number of cells drawn per group in each block,
#'                  picked at random (`seed`), default is 200. NULL draws all
#'                  cells. Only for `columns = "cells"`
#'
#' @param scale Whether to z-score each feature, default is TRUE: across all
#'              plotted cells (before downsampling) for `columns = "cells"`,
#'              across the columns for `columns = "groups"`
#'
#' @param min.cutoff,max.cutoff Ends of the colour scale: a number, a
#'                              quantile of the plotted values such as "q5"
#'                              / "q95", or NA (default). For z-scores, NA
#'                              gives +-min(largest |z|, 2.5); for values
#'                              that are all positive, 0 and the largest
#'                              value; for other values, +-largest |value|.
#'                              Values beyond the ends get the colour of the
#'                              end
#'
#' @param group.marker How the groups are marked above the columns: "dot"
#'                     (default), "strip" or "none"
#'
#' @param group.cols Colours of the groups: NULL (default) for the colours
#'                   \code{\link{DimPlot_bioToolkit}} gives them, a
#'                   qualitative palette name, or colours named by group
#'
#' @param annotations Further discrete metadata columns shown as colour strips
#'                    above the columns (e.g. "sample"), each with its legend.
#'                    With `columns = "groups"`, each must have one value per
#'                    column. Cells in a group are ordered by the first one
#'
#' @param annotation.cols Colours of `annotations`: a list named by column,
#'                        each a qualitative palette name or colours named by
#'                        value. Columns not in the list get "tol-muted",
#'                        "okabe-ito", "tol-light", ... in turn (or the default
#'                        group colours when they have too many values)
#'
#' @param cols Colours of the heatmap. If NULL (default), for z-scores or
#'             signed values: "teal-yellow" (teal - black - yellow) with
#'             `columns = "cells"`, "blue-magenta" (light centre) with
#'             `columns = "groups"`; "magma" for unscaled non-negative values. Otherwise a palette
#'             name (see \code{\link{show_palettes}}) or colours from low to
#'             high, e.g. `c(low, mid, high)`
#'
#' @param label.features Which feature names are drawn: NULL (default) all
#'                       of them if a row is tall enough for the text (at
#'                       least 0.125 in), else none; TRUE or FALSE for all or
#'                       none; or a few features, named with leader lines
#'
#' @param label.side Side of the feature names: "right" (default) or "left".
#'                   Names of blocks of features are on the other side
#'
#' @param italic.genes Whether gene names are set in italics, default is TRUE
#'                     (genomic regions such as peaks never are)
#'
#' @param panel.width Width of the heatmap in inches for `columns = "cells"`
#'                    (without the gaps between blocks), default is 6
#'
#' @param cell.width Width of one column in inches for `columns = "groups"`,
#'                   default is 0.18
#'
#' @param cell.height Height of one row in inches, default is 0.13, which
#'                    fits the names. Used while the heatmap is at most 6 in
#'                    tall; beyond that, see `panel.height`
#'
#' @param panel.height Height of the heatmap in inches (without the gaps
#'                     between blocks of features). If NULL (default), rows
#'                     of `cell.height`, up to 6 in in all
#'
#' @param split.gap Gap between blocks of columns in inches, default is 0.04
#'
#' @param seed Random seed of the downsampling and cell order, default is 42
#'
#' @param legend.position Position of the legends: "right" (default),
#'                        "bottom", "top", "left" or "none"
#'
#' @param bg Background colour, default is "transparent"
#'
#' @details
#' Cells without a value of `group.by` or `split.by` are left out, with a
#' message. A feature with the same value in every cell has no z-score and is
#' drawn at 0, with a message.
#'
#' The source data (\code{\link{export_source_data}}) hold one row per column
#' and feature: the `column` (cell, or group mean), `block`, `group`,
#' `split`, the `annotations`, `feature_group`, `feature`, `value` (as read,
#' or the group mean), `value_scaled` (when scaled), `n_cells` (group means)
#' and the `assay`, `layer` and `average` used.
#'
#' @return A patchwork (body, header rows and labels)
#'
#' @examples
#' \dontrun{
#' markers <- top_markers(Seurat::FindAllMarkers(obj, only.pos = TRUE), n = 5)
#'
#' # Cells, z-scored per gene, at most 200 per cell type
#' DoHeatmap_bioToolkit(obj, markers, group.by = "cell_type", assay = "SCT")
#'
#' # Group means: conditions (dots) compared within each cell type (blocks)
#' p <- DoHeatmap_bioToolkit(obj, markers, group.by = "condition",
#'                           split.by = "cell_type", columns = "groups")
#' save_plot(p, "heatmap.pdf")
#'
#' # Many genes: names hidden, a few named with leader lines
#' DoHeatmap_bioToolkit(obj, top_markers(all_markers, n = 30), group.by = "cell_type",
#'                      label.features = c("Col1a1", "Cdh5", "Myh6"))
#' }
#'
#' @seealso \code{\link{top_markers}}, \code{\link{DotPlot_bioToolkit}},
#'          \code{\link{save_plot}}, \code{\link{export_source_data}}
#'
#' @import ggplot2
#' @export
#'
DoHeatmap_bioToolkit <- function(
    object,
    features,
    group.by = NULL,
    split.by = NULL,
    cells = NULL,
    assay = NULL,
    slot = "data",
    columns = c("cells", "groups"),
    average = c("linear", "log"),
    max.cells = 200,
    scale = TRUE,
    min.cutoff = NA,
    max.cutoff = NA,
    group.marker = c("dot", "strip", "none"),
    group.cols = NULL,
    annotations = NULL,
    annotation.cols = NULL,
    cols = NULL,
    label.features = NULL,
    label.side = c("right", "left"),
    italic.genes = TRUE,
    panel.width = 6,
    cell.width = 0.18,
    cell.height = 0.13,
    panel.height = NULL,
    split.gap = 0.04,
    seed = 42,
    legend.position = c("right", "bottom", "top", "left", "none"),
    bg = "transparent"
) {
    check_installed("Seurat", reason = "to fetch expression")
    columns <- match.arg(columns)
    average <- match.arg(average)
    group.marker <- match.arg(group.marker)
    label.side <- match.arg(label.side)
    legend.position <- match.arg(legend.position)
    if (!is.null(max.cells) && !(is.numeric(max.cells) && length(max.cells) == 1 && max.cells >= 1)) {
        stop("`max.cells` must be NULL or a single number >= 1", call. = FALSE)
    }
    sz <- theme_sizes()

    # Cells, their groups, blocks and annotations ---------------------------------
    hd <- heatmap_cells(object, group.by, split.by, annotations, cells)

    # Values: cells x features ------------------------------------------------------
    ft <- dot_feature_table(features)
    fv <- fetch_plot_values(object, unique(ft$feature), hd$cells$cell, assay, slot)
    ft <- ft[ft$feature %in% colnames(fv$data), , drop = FALSE]
    expr <- as.matrix(fv$data)
    is_gene <- fv$is_gene
    if (slot == "scale.data" && scale) {
        message("Layer 'scale.data' is drawn as stored (already scaled)")
        scale <- FALSE
    }
    sv <- heatmap_values(expr, hd, columns, average, scale, is_gene, log.data = slot == "data", max.cells, seed)

    # Colours and limits --------------------------------------------------------------
    fill <- sv$fill[, ft$feature, drop = FALSE]
    signed <- scale || slot == "scale.data" || any(fill < 0, na.rm = TRUE)
    limits <- heatmap_fill_limits(as.vector(fill), scale, signed, min.cutoff, max.cutoff)
    fill_cols <- resolve_continuous_palette(cols, diverging = FALSE,
                                            default = if (!signed) "magma" else if (columns == "cells") "teal-yellow" else "blue-magenta")
    any_gene <- any(is_gene[ft$feature])
    fill_title <- paste0(if (any_gene) "Expression" else "Value",
                         if (scale) "\n(z-score)" else if (slot == "scale.data") "\n(scaled)")

    # Layout: columns and rows -----------------------------------------------------------
    cl <- heatmap_column_layout(sv$columns, columns, panel.width, cell.width)
    rl <- heatmap_row_layout(ft, cell.height, panel.height)
    show_names <- if (is.null(label.features)) {
        rl$row.height >= sz$small / 72
    } else if (is.logical(label.features)) {
        isTRUE(label.features)
    } else {
        FALSE
    }
    if (is.null(label.features) && !show_names) {
        message("Rows are ", signif(rl$row.height, 2), " in tall, too thin to name each feature: names ",
                "are hidden (`label.features` names a few, `cell.height` / `panel.height` enlarge rows)")
    }
    mark <- if (is.character(label.features)) {
        missing <- setdiff(label.features, ft$feature)
        if (length(missing)) warning("`label.features` not plotted: ", paste(missing, collapse = ", "),
                                     call. = FALSE)
        intersect(label.features, ft$feature)
    } else {
        character(0)
    }
    italic <- italic.genes & is_gene

    # Plots ----------------------------------------------------------------------------------
    long <- heatmap_long(sv, ft, cl, rl)
    body <- heatmap_body(long, cl, rl, fill_cols, limits, fill_title, show_names, italic, label.side,
                         split.gap, legend.position, sz, bg)
    group_all <- hd$group_levels
    group_cols <- group_colours(group.cols, group_all, arg = "group.cols")$cols
    header <- list()
    if (group.marker != "none") {
        header$group <- heatmap_marker_row(cl, group_cols, group.marker, columns,
                                           legend = !is.null(split.by),
                                           title = if (is.null(group.by)) "Identity" else group.by,
                                           split.gap, sz, bg)
    }
    for (i in seq_along(annotations)) {
        a <- annotations[i]
        header[[a]] <- heatmap_annotation_row(cl, a, annotation_colours(annotation.cols, a, hd$cells[[a]], i),
                                              split.gap, sz, bg)
    }
    # Block names on the top row
    rows <- c(header, list(body = body))
    rows[[1]] <- rows[[1]] + heatmap_block_labels(cl, sz)
    for (i in seq_along(rows)[-1]) rows[[i]] <- rows[[i]] + theme(strip.text.x = element_blank())

    side <- if (length(mark)) heatmap_mark_labels(rl, ft, mark, italic, label.side, sz, bg)
    p <- heatmap_assemble(rows, side, label.side, legend.position, bg)

    # Source data: on the body
    src <- long$source
    if (!is.null(fv$assay)) {
        src$assay <- fv$assay$assay
        src$layer <- slot
    }
    if (columns == "groups") src$average <- if (slot == "data") average else "arithmetic"
    attr(p, "source_data") <- src
    p
}


#' Cells of a heatmap with their group, block and annotations
#'
#' @param object A Seurat object
#'
#' @param group.by,split.by,annotations,cells As in
#'        \code{\link{DoHeatmap_bioToolkit}}
#'
#' @return A list with `cells` (a data frame: `cell`, `group`, `split` or NA,
#'         `block`, and one column per annotation) and `group_levels` (a
#'         factor with all levels of `group.by`, for the colours)
#'
#' @keywords internal
heatmap_cells <- function(object,
                          group.by,
                          split.by,
                          annotations,
                          cells) {
    cells <- resolve_cells(cells, colnames(object))
    if (length(cells) == 0) stop("No cells to plot", call. = FALSE)
    group_all <- fetch_group(object, group.by, colnames(object))
    group <- group_all[match(cells, colnames(object))]
    split <- NULL
    if (!is.null(split.by)) {
        if (length(split.by) != 1) stop("`split.by` must be a single metadata column", call. = FALSE)
        split <- fetch_group(object, split.by, cells)
    }
    ann <- lapply(annotations, function(a) {
        v <- fetch_meta(object, a, cells)
        if (is.numeric(v) && !is.factor(v)) {
            stop("Annotation '", a, "' is numeric; annotations must be discrete (e.g. sample)", call. = FALSE)
        }
        if (is.factor(v)) v else factor(v)
    })
    names(ann) <- annotations
    keep <- !is.na(group) & (if (is.null(split)) TRUE else !is.na(split))
    if (!all(keep)) {
        message(sum(!keep), " cell(s) without a value of `",
                paste(c(if (is.null(group.by)) "identity" else group.by, split.by), collapse = "` or `"),
                "` are left out")
    }
    if (!any(keep)) stop("No cells to plot", call. = FALSE)
    out <- data.frame(cell = cells[keep], stringsAsFactors = FALSE)
    out$group <- droplevels(group[keep])
    out$split <- if (is.null(split)) factor(NA) else droplevels(split[keep])
    # Blocks: the splits, or the groups themselves
    out$block <- if (is.null(split)) out$group else out$split
    for (a in annotations) out[[a]] <- droplevels(ann[[a]][keep])
    list(cells = out, group_levels = group_all, annotations = annotations)
}


#' Values of a heatmap: z-scored cells (downsampled) or group means
#'
#' @param expr Matrix of cells x features, rows in the order of `hd$cells`
#'
#' @param hd Result of \code{heatmap_cells()}
#'
#' @param columns "cells" or "groups"
#'
#' @param average "linear" or "log"
#'
#' @param scale Whether to z-score each feature
#'
#' @param is_gene Logical, named by feature
#'
#' @param log.data Whether genes hold log1p values (the "data" layer)
#'
#' @param max.cells Largest number of cells per group and block, or NULL
#'
#' @param seed Random seed
#'
#' @return A list with `columns` (a data frame, one row per column: `column`,
#'         `block`, `group`, `split`, annotations and `n_cells`), `value`
#'         and `fill` (matrices of columns x features: values as read or
#'         averaged, and as coloured)
#'
#' @keywords internal
heatmap_values <- function(expr,
                           hd,
                           columns,
                           average,
                           scale,
                           is_gene,
                           log.data,
                           max.cells,
                           seed) {
    cd <- hd$cells
    rownames(expr) <- cd$cell
    if (scale) {
        sds <- apply(expr, 2, stats::sd)
        flat <- names(sds)[!is.finite(sds) | sds == 0]
        if (length(flat)) message("Same value in every cell, drawn at 0: ", paste(flat, collapse = ", "))
    }
    unit <- interaction(cd$block, cd$group, drop = TRUE, lex.order = TRUE, sep = " | ")

    if (columns == "cells") {
        # z-scores from all plotted cells, before downsampling
        fill <- if (scale) scale_columns(expr) else expr
        # Downsample per group and block, then order: blocks, groups, then by
        # the first annotation (random order within it)
        keep <- with_seed(seed, {
            idx <- split(seq_len(nrow(cd)), unit, drop = TRUE)
            unlist(lapply(idx, function(i) {
                if (!is.null(max.cells) && length(i) > max.cells) i <- i[sample.int(length(i), max.cells)]
                i[sample.int(length(i))]
            }), use.names = FALSE)
        })
        first_ann <- if (length(hd$annotations)) cd[[hd$annotations[1]]][keep] else rep(1, length(keep))
        keep <- keep[order(cd$block[keep], cd$group[keep], first_ann)]
        col_df <- cd[keep, , drop = FALSE]
        col_df$column <- col_df$cell
        col_df$n_cells <- 1L
        rownames(col_df) <- NULL
        return(list(columns = col_df, value = expr[keep, , drop = FALSE], fill = fill[keep, , drop = FALSE]))
    }

    # Group means: one column per group in each block
    linear <- log.data & average == "linear" & is_gene[colnames(expr)]
    x <- expr
    x[, linear] <- expm1(x[, linear])
    n <- as.vector(table(unit))
    mean <- rowsum(x, unit, reorder = TRUE) / n
    mean[, linear] <- log1p(mean[, linear])
    lv <- levels(unit)
    first <- match(lv, as.character(unit))
    col_df <- cd[first, c("block", "group", "split", hd$annotations), drop = FALSE]
    # Annotations must have one value per column
    for (a in hd$annotations) {
        per <- tapply(as.character(cd[[a]]), unit, function(v) length(unique(v)))
        if (any(per > 1, na.rm = TRUE)) {
            stop("Annotation '", a, "' has several values within a column (group in a block); with ",
                 "columns = \"groups\", annotations must be constant within each (e.g. condition ",
                 "per sample)", call. = FALSE)
        }
    }
    col_df$column <- lv
    col_df$n_cells <- n
    rownames(col_df) <- NULL
    rownames(mean) <- lv
    fill <- if (scale) {
        if (nrow(mean) < 2) {
            message("Only one column: values are not scaled")
            mean
        } else {
            scale_columns(mean)
        }
    } else {
        mean
    }
    list(columns = col_df, value = mean, fill = fill)
}


#' Colour limits of a heatmap
#'
#' @param v Values coloured
#'
#' @param scale Whether they are z-scores
#'
#' @param signed Whether the scale is centred on 0
#'
#' @param min.cutoff,max.cutoff Cutoffs (number, "q##" or NA)
#'
#' @return Two numbers
#'
#' @keywords internal
heatmap_fill_limits <- function(v,
                                scale,
                                signed,
                                min.cutoff,
                                max.cutoff) {
    v <- v[is.finite(v)]
    if (length(v) == 0) return(c(-1, 1))
    if (scale || !signed) return(dot_fill_limits(v, scale, min.cutoff, max.cutoff))
    # Signed values as they are: symmetric around 0
    m <- max(abs(v))
    if (m == 0) m <- 1
    lim <- c(cutoff_value(min.cutoff, v, function(x) -m), cutoff_value(max.cutoff, v, function(x) m))
    if (lim[2] <= lim[1]) lim[2] <- lim[1] + 1
    lim
}


#' Positions and sizes of the columns of a heatmap
#'
#' Columns sit at integer positions within each block (facet). For cells,
#' groups within a block are set apart by empty columns about 0.015 in wide.
#'
#' @param col_df Column table from \code{heatmap_values()}
#'
#' @param columns "cells" or "groups"
#'
#' @param panel.width Width of all columns of cells (inches)
#'
#' @param cell.width Width of a column of group means (inches)
#'
#' @return A list with `columns` (`col_df` with `x`), `widths` (inches per
#'         block), `width` (one column, inches), `limits` (x range per block)
#'         and `blocks` (block levels)
#'
#' @keywords internal
heatmap_column_layout <- function(col_df,
                                  columns,
                                  panel.width,
                                  cell.width) {
    blocks <- levels(droplevels(col_df$block))
    col_df$block <- factor(as.character(col_df$block), levels = blocks)
    n_groups <- vapply(blocks, function(b) length(unique(col_df$group[col_df$block == b])), integer(1))
    if (columns == "cells") {
        # Gap columns between groups within a block, so that 0.015 in gaps
        # stay a whole number of columns
        n_inner <- sum(n_groups - 1)
        width <- panel.width / nrow(col_df)
        gap <- if (n_inner > 0) max(1L, as.integer(round(0.015 / width))) else 0L
        width <- panel.width / (nrow(col_df) + n_inner * gap)
    } else {
        width <- cell.width
        gap <- 0L
    }
    col_df$x <- NA_real_
    for (b in blocks) {
        i <- which(col_df$block == b)
        g <- as.integer(droplevels(col_df$group[i]))
        # Position: running index plus the gaps before the column's group
        col_df$x[i] <- seq_along(i) + (g - 1L) * gap
    }
    span <- vapply(blocks, function(b) max(col_df$x[col_df$block == b]), numeric(1))
    list(
        columns = col_df,
        widths = span * width,
        width = width,
        limits = lapply(stats::setNames(span, blocks), function(s) c(0.5, s + 0.5)),
        blocks = blocks
    )
}


#' Positions and sizes of the rows of a heatmap
#'
#' @param ft Feature table from \code{dot_feature_table()}
#'
#' @param cell.height Height of one row (inches)
#'
#' @param panel.height Height of all rows (inches), or NULL: rows of
#'                     `cell.height`, up to 6 in
#'
#' @return A list with `rows` (`ft` with `y`, 1 at the bottom, and `block`),
#'         `heights` (inches per block, top to bottom), `row.height` and
#'         `blocks`
#'
#' @keywords internal
heatmap_row_layout <- function(ft,
                               cell.height,
                               panel.height) {
    n <- nrow(ft)
    if (is.null(panel.height)) panel.height <- min(n * cell.height, 6)
    row_h <- panel.height / n
    ft$block <- if (is.null(ft$feature_group)) factor("") else ft$feature_group
    # First feature at the top
    ft$y <- rev(seq_len(n))
    blocks <- levels(droplevels(ft$block))
    heights <- vapply(blocks, function(b) sum(ft$block == b), integer(1)) * row_h
    list(rows = ft, heights = heights, row.height = row_h, blocks = blocks,
         grouped = !is.null(ft$feature_group))
}


#' Long table of heatmap tiles, and the source data
#'
#' @param sv Result of \code{heatmap_values()}
#'
#' @param ft Feature table (rows of the heatmap)
#'
#' @param cl Result of \code{heatmap_column_layout()}
#'
#' @param rl Result of \code{heatmap_row_layout()}
#'
#' @return A list with `tiles` (x, y, fill, block, row_block) and `source`
#'
#' @keywords internal
heatmap_long <- function(sv,
                         ft,
                         cl,
                         rl) {
    cols <- cl$columns
    rows <- rl$rows
    idx <- expand.grid(i = seq_len(nrow(cols)), j = seq_len(nrow(rows)))
    fill <- sv$fill[cbind(idx$i, match(rows$feature[idx$j], colnames(sv$fill)))]
    value <- sv$value[cbind(idx$i, match(rows$feature[idx$j], colnames(sv$value)))]
    tiles <- data.frame(
        x = cols$x[idx$i],
        y = rows$y[idx$j],
        fill = fill,
        block = cols$block[idx$i],
        row_block = factor(as.character(rows$block[idx$j]), levels = rl$blocks)
    )
    ann <- setdiff(names(cols), c("cell", "column", "block", "group", "split", "n_cells", "x"))
    src <- data.frame(column = cols$column[idx$i], block = cols$block[idx$i], group = cols$group[idx$i],
                      stringsAsFactors = FALSE)
    if (!all(is.na(cols$split))) src$split <- cols$split[idx$i]
    for (a in ann) src[[a]] <- cols[[a]][idx$i]
    if (rl$grouped) src$feature_group <- rows$feature_group[idx$j]
    src$feature <- rows$feature[idx$j]
    src$value <- value
    if (!identical(sv$fill, sv$value)) src$value_scaled <- fill
    if (any(cols$n_cells > 1)) src$n_cells <- cols$n_cells[idx$i]
    list(tiles = tiles, source = src)
}


#' Facets shared by every row of a heatmap
#'
#' @param rows Whether to facet by blocks of rows too
#'
#' @param label.side Side of the feature names; block names are on the other
#'
#' @return A facet
#'
#' @keywords internal
heatmap_facet <- function(rows,
                          label.side = "right") {
    if (rows) {
        facet_grid(rows = vars(.data$row_block), cols = vars(.data$block), scales = "free",
                   space = "free", switch = if (label.side == "right") "y" else NULL)
    } else {
        facet_grid(cols = vars(.data$block), scales = "free_x", space = "free_x")
    }
}


#' The x scale and theme shared by every row of a heatmap
#'
#' @param cl Result of \code{heatmap_column_layout()}
#'
#' @param height Panel heights (inches)
#'
#' @param split.gap Gap between blocks of columns (inches)
#'
#' @param bg Background colour
#'
#' @return A list of ggplot2 components
#'
#' @keywords internal
heatmap_frame <- function(cl,
                          height,
                          split.gap,
                          bg) {
    list(
        scale_x_continuous(expand = expansion(0), breaks = NULL),
        theme_bioToolkit("grid", bg = bg),
        panel_size_theme(unname(cl$widths), unname(height), panel.spacing = split.gap),
        theme(
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            panel.grid.major = element_blank(),
            strip.placement = "outside",
            strip.clip = "off",
            plot.margin = margin(1, 2, 1, 2)
        )
    )
}


#' Body of a heatmap: the tiles, feature names and blocks of rows
#'
#' @keywords internal
heatmap_body <- function(long,
                         cl,
                         rl,
                         fill_cols,
                         limits,
                         fill_title,
                         show_names,
                         italic,
                         label.side,
                         split.gap,
                         legend.position,
                         sz,
                         bg) {
    tiles <- long$tiles
    rows <- rl$rows
    # Range of each block of rows: blank points at its ends
    ends <- do.call(rbind, lapply(rl$blocks, function(b) {
        y <- rows$y[rows$block == b]
        data.frame(y = c(min(y) - 0.5, max(y) + 0.5), row_block = factor(b, levels = rl$blocks))
    }))
    ends <- merge(ends, data.frame(block = factor(cl$blocks, levels = cl$blocks)))
    ends$x <- vapply(as.character(ends$block), function(b) cl$limits[[b]][1], numeric(1))
    ends2 <- ends
    ends2$x <- vapply(as.character(ends$block), function(b) cl$limits[[b]][2], numeric(1))
    labels <- stats::setNames(rows$feature, rows$y)
    vertical <- !(legend.position %in% c("top", "bottom"))
    # Gap columns between groups hold transparent tiles, so that the raster
    # grid of each panel stays regular
    gaps <- do.call(rbind, lapply(cl$blocks, function(b) {
        used <- cl$columns$x[cl$columns$block == b]
        free <- setdiff(seq_len(max(used)), used)
        if (length(free) == 0) return(NULL)
        g <- expand.grid(x = free, j = seq_len(nrow(rows)))
        data.frame(x = g$x, y = rows$y[g$j], fill = NA_real_, block = factor(b, levels = cl$blocks),
                   row_block = factor(as.character(rows$block[g$j]), levels = rl$blocks))
    }))
    if (!is.null(gaps)) tiles <- rbind(tiles, gaps)
    p <- ggplot(tiles, aes(x = .data$x, y = .data$y, fill = .data$fill)) +
        geom_blank(data = rbind(ends, ends2), inherit.aes = FALSE, aes(x = .data$x, y = .data$y)) +
        geom_raster() +
        scale_fill_gradientn(
            colours = fill_cols, limits = limits, oob = scales::oob_squish,
            rescaler = function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
                if (limits[1] < 0 && limits[2] > 0) scales::rescale_mid(x, to, from, mid = 0)
                else scales::rescale(x, to, from)
            },
            na.value = "transparent", name = fill_title,
            guide = guide_colourbar(
                order = 1,
                direction = if (vertical) "vertical" else "horizontal",
                theme = theme(
                    legend.key.height = unit(if (vertical) 0.9 else 0.12, "in"),
                    legend.key.width = unit(if (vertical) 0.12 else 0.9, "in"),
                    legend.frame = element_rect(colour = "black", linewidth = 0.3),
                    legend.ticks = element_line(colour = "black", linewidth = 0.3),
                    legend.ticks.length = unit(0.03, "in")
                )
            )
        ) +
        scale_y_continuous(
            expand = expansion(0),
            breaks = if (show_names) rows$y else NULL,
            labels = function(b) {
                f <- unname(labels[as.character(b)])
                italic_labels(f, unname(italic[f]))
            },
            position = label.side
        ) +
        heatmap_facet(TRUE, label.side) +
        heatmap_frame(cl, rl$heights, split.gap, bg) +
        theme(
            axis.text.y = element_text(size = sz$small, margin = margin(l = 2, r = 2)),
            panel.spacing.y = unit(0.04, "in"),
            strip.text.y.left = element_text(size = sz$small, face = "bold", angle = 0, hjust = 1,
                                             margin = margin(r = 3)),
            strip.text.y.right = element_text(size = sz$small, face = "bold", angle = 0, hjust = 0,
                                              margin = margin(l = 3)),
            legend.title = element_text(size = sz$small, margin = margin(b = 0.15, unit = "cm")),
            legend.title.position = "top"
        )
    if (!rl$grouped) p <- p + theme(strip.text.y = element_blank())
    p
}


#' Header row of a heatmap marking the groups with their colours
#'
#' @keywords internal
heatmap_marker_row <- function(cl,
                               group_cols,
                               style,
                               columns,
                               legend,
                               title,
                               split.gap,
                               sz,
                               bg) {
    cols <- cl$columns
    # One mark per group in each block, over the group's columns
    key <- paste(cols$block, cols$group, sep = "\r")
    span <- cols[!duplicated(key), c("block", "group"), drop = FALSE]
    k <- paste(span$block, span$group, sep = "\r")
    span$xmin <- unname(tapply(cols$x, key, min)[k]) - 0.5
    span$xmax <- unname(tapply(cols$x, key, max)[k]) + 0.5
    span$x <- (span$xmin + span$xmax) / 2
    span$group <- factor(as.character(span$group), levels = names(group_cols))
    span$y <- 0
    height <- if (style == "dot") 0.16 else 0.09
    # Dot diameter: 70% of a column of group means, else a fixed 0.11 in
    dot_in <- if (columns == "groups") min(0.7 * cl$width, 0.11) else 0.11
    # A point of size s is about 0.75 * s * .pt points across (as in the dot grid)
    dot_size <- dot_in * 72 / (0.75 * ggplot2::.pt)
    p <- ggplot() + heatmap_blank_blocks(cl)
    p <- p + if (style == "dot") {
        geom_point(data = span, aes(x = .data$x, y = .data$y, colour = .data$group), shape = 16,
                   size = dot_size)
    } else {
        geom_rect(data = span, aes(xmin = .data$xmin, xmax = .data$xmax, ymin = -0.5, ymax = 0.5,
                                   fill = .data$group), colour = NA)
    }
    guide <- if (legend) guide_legend(order = 2, override.aes = list(size = 3)) else "none"
    p <- p +
        (if (style == "dot") {
            scale_colour_manual(values = group_cols, drop = TRUE, name = title, guide = guide)
        } else {
            scale_fill_manual(values = group_cols, drop = TRUE, name = title,
                              guide = if (legend) guide_legend(order = 2) else "none")
        }) +
        scale_y_continuous(limits = c(-0.5, 0.5), expand = expansion(0), breaks = NULL) +
        heatmap_facet(FALSE) +
        heatmap_frame(cl, height, split.gap, bg) +
        theme(legend.title = element_text(size = sz$small), legend.key.size = unit(0.16, "in"))
    p
}


#' Header row of a heatmap: a colour strip of one annotation
#'
#' @keywords internal
heatmap_annotation_row <- function(cl,
                                   annotation,
                                   colours,
                                   split.gap,
                                   sz,
                                   bg) {
    cols <- cl$columns
    d <- data.frame(x = cols$x, y = 0, block = cols$block, value = cols[[annotation]])
    ggplot() +
        heatmap_blank_blocks(cl) +
        geom_tile(data = d, aes(x = .data$x, y = .data$y, fill = .data$value), width = 1, height = 1,
                  colour = NA) +
        scale_fill_manual(values = colours, drop = TRUE, name = annotation,
                          guide = guide_legend(order = 3)) +
        scale_y_continuous(limits = c(-0.5, 0.5), expand = expansion(0), breaks = NULL) +
        heatmap_facet(FALSE) +
        heatmap_frame(cl, 0.09, split.gap, bg) +
        theme(legend.title = element_text(size = sz$small), legend.key.size = unit(0.16, "in"))
}


#' Blank points at the x range of every block, so that header rows match the
#' body
#'
#' @keywords internal
heatmap_blank_blocks <- function(cl) {
    d <- data.frame(
        block = factor(rep(cl$blocks, each = 2), levels = cl$blocks),
        x = unlist(cl$limits, use.names = FALSE),
        y = 0
    )
    geom_blank(data = d, aes(x = .data$x, y = .data$y))
}


#' Names of the blocks of columns, as the facet strips of the top row
#'
#' Turned upright when a name is wider than its block (estimated from the
#' number of characters, as in the dot grid).
#'
#' @keywords internal
heatmap_block_labels <- function(cl,
                                 sz) {
    label_w <- nchar(cl$blocks) * 0.6 * sz$small / 72
    upright <- any(label_w > cl$widths + 0.04)
    theme(strip.text.x.top = element_text(size = sz$small, face = "bold",
                                          angle = if (upright) 90 else 0,
                                          hjust = if (upright) 0 else 0.5,
                                          vjust = 0.5,
                                          margin = margin(b = 3)))
}


#' Colours of an annotation strip
#'
#' @param annotation.cols List of colours named by annotation, or NULL
#'
#' @param a Annotation name
#'
#' @param values Its values (factor)
#'
#' @param i Index of the annotation, choosing its default palette
#'
#' @return Colours named by value
#'
#' @keywords internal
annotation_colours <- function(annotation.cols,
                               a,
                               values,
                               i) {
    values <- droplevels(values)
    given <- if (!is.null(annotation.cols)) annotation.cols[[a]]
    if (!is.null(given)) return(group_colours(given, values, arg = paste0("annotation.cols$", a))$cols)
    defaults <- c("tol-muted", "okabe-ito", "tol-light", "tol-bright", "tol-vibrant")
    pal <- defaults[(i - 1) %% length(defaults) + 1]
    if (nlevels(values) > length(bt_palette(pal))) pal <- NULL
    group_colours(pal, values)$cols
}


#' A few feature names beside the heatmap, with leader lines
#'
#' Names are spread apart so that they do not overlap (at least one text
#' line apart, within each block of rows), and joined to their rows by
#' short lines, as ComplexHeatmap's anno_mark().
#'
#' @keywords internal
heatmap_mark_labels <- function(rl,
                                ft,
                                mark,
                                italic,
                                label.side,
                                sz,
                                bg) {
    rows <- rl$rows
    line_rows <- (sz$small / 72 * 1.15) / rl$row.height
    marks <- rows[rows$feature %in% mark, , drop = FALSE]
    marks$row_block <- factor(as.character(marks$block), levels = rl$blocks)
    marks$label_y <- NA_real_
    for (b in unique(as.character(marks$row_block))) {
        i <- which(marks$row_block == b)
        lim <- range(rows$y[rows$block == b]) + c(-0.5, 0.5)
        marks$label_y[i] <- spread_positions(marks$y[i], line_rows, lim)
    }
    # Labels as plotmath strings, parsed by geom_text()
    quoted <- vapply(marks$feature, deparse, character(1))
    marks$label <- ifelse(italic[marks$feature], paste0("italic(", quoted, ")"), quoted)
    # Panel: leader lines (0.25 in), then the names; its width holds the
    # longest name (about 0.55 em per character), so that names never reach
    # into the legends
    room <- max(nchar(marks$feature)) * 0.55 * sz$small / 72
    width <- 0.25 + room + 0.05
    right <- label.side == "right"
    seg <- data.frame(row_block = marks$row_block, y = marks$y, label_y = marks$label_y)
    ends <- do.call(rbind, lapply(rl$blocks, function(b) {
        y <- rows$y[rows$block == b]
        data.frame(y = c(min(y) - 0.5, max(y) + 0.5), row_block = factor(b, levels = rl$blocks), x = 0.1)
    }))
    # x in inches from the heatmap edge
    x0 <- if (right) 0.02 else width - 0.02
    x1 <- if (right) 0.08 else width - 0.08
    x2 <- if (right) 0.2 else width - 0.2
    ggplot() +
        geom_blank(data = ends, aes(x = .data$x, y = .data$y)) +
        geom_segment(data = seg, aes(x = x0, xend = x1, y = .data$y, yend = .data$y),
                     linewidth = pt_to_linewidth(sz$line), colour = "grey30") +
        geom_segment(data = seg, aes(x = x1, xend = x2, y = .data$y, yend = .data$label_y),
                     linewidth = pt_to_linewidth(sz$line), colour = "grey30") +
        geom_text(data = marks, aes(x = x2 + if (right) 0.03 else -0.03, y = .data$label_y,
                                    label = .data$label),
                  parse = TRUE, hjust = if (right) 0 else 1, size = pt_to_size(sz$small)) +
        scale_x_continuous(limits = c(0, width), expand = expansion(0)) +
        scale_y_continuous(expand = expansion(0)) +
        facet_grid(rows = vars(.data$row_block), scales = "free_y", space = "free_y") +
        coord_cartesian(clip = "off") +
        theme_void() +
        panel_size_theme(width, unname(rl$heights), panel.spacing = 0.04) +
        theme(strip.text = element_blank(),
              plot.margin = margin(1, 2, 1, 2),
              plot.background = theme_bg(bg)$plot.background)
}


#' Spread positions apart to a minimum distance, within limits
#'
#' @param y Positions
#'
#' @param gap Smallest distance between neighbours
#'
#' @param lim Range the positions are kept in, if they fit
#'
#' @return The new positions, in the order of `y`
#'
#' @keywords internal
spread_positions <- function(y,
                             gap,
                             lim) {
    o <- order(y)
    s <- y[o]
    n <- length(s)
    if (n > 1) {
        # Push upwards, then shift the whole set back down into the limits
        # and push downwards, so that labels stay near their rows
        for (i in 2:n) s[i] <- max(s[i], s[i - 1] + gap)
        over <- s[n] - (lim[2] - gap / 2)
        if (over > 0) s <- s - over
        for (i in (n - 1):1) s[i] <- min(s[i], s[i + 1] - gap)
    }
    out <- numeric(n)
    out[o] <- s
    out
}


#' Join the rows of a heatmap and the labels with leader lines
#'
#' @keywords internal
heatmap_assemble <- function(rows,
                             side,
                             label.side,
                             legend.position,
                             bg) {
    n <- length(rows)
    # Spacers are drawn on white by default; keep the background
    spacer <- function() patchwork::plot_spacer() + theme_bg(bg)
    if (is.null(side)) {
        plots <- rows
        layout <- patchwork::plot_layout(ncol = 1, guides = "collect")
    } else {
        # Header rows over the body only; the labels beside the body
        if (label.side == "right") {
            plots <- c(unlist(lapply(rows[-n], function(r) list(r, spacer())), recursive = FALSE),
                       list(rows[[n]], side))
        } else {
            plots <- c(unlist(lapply(rows[-n], function(r) list(spacer(), r)), recursive = FALSE),
                       list(side, rows[[n]]))
        }
        layout <- patchwork::plot_layout(ncol = 2, guides = "collect")
    }
    patchwork::wrap_plots(plots) + layout +
        patchwork::plot_annotation(theme = theme_bg(bg) + theme(legend.position = legend.position))
}
