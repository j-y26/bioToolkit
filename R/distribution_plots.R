# Per-cell distributions: violin and ridge plots (roadmap item S1)
#
# Design
# ------
# `VlnPlot_bioToolkit()` and `RidgePlot_bioToolkit()` replace
# Seurat::VlnPlot() / RidgePlot() and the user's own QC violin and ridge
# plots (`qc_violin()`, `qc_ridge()`), whose look is kept as the default:
# violins without outline, scaled to the same width; a white box of the
# inter-quartile range and median, without whiskers or outliers; a random
# sample of the cells as rasterized points; ridges with a white outline.
#
# * Any feature: genes (read from one assay and layer, reported in a message;
#   R/seurat_data.R), metadata columns (QC metrics, module scores) and
#   embedding dimensions, mixed in one call. `object` may also be a data
#   frame (e.g. a QC table before a Seurat object exists).
# * One plot per feature, combined with patchwork, so that each feature gets
#   its own value axis: `log` and `percent` take TRUE or the names of the
#   features they apply to, since one call often mixes counts (log axis),
#   fractions (% axis) and expression (linear axis).
# * `y.max` (Seurat's name) zooms the value axis to a number or a quantile
#   ("q99") of the values: densities and boxes are still computed from all
#   cells (coord_cartesian()), unlike scale limits, which drop cells. The
#   panels are not clipped (the split.by bars sit outside them), so violins
#   and boxes are cut at the top themselves (`GeomViolinTop`,
#   `GeomBoxplotTop`) and points above it are not drawn (2026-10-14; they
#   used to run past the axis).
# * Genes are sparse: most cells of a group may have 0, and the violin is
#   then a spike at 0. `zeros = "exclude"` draws the cells above 0 only;
#   "annotate" also prints the % of cells above 0 for each group.
# * `group.by` / `split.by` (package-wide rule, 2026-10-13): `group.by` gives
#   the violins, compared with each other and filled with the group colours
#   (as in DimPlot_bioToolkit()); `split.by` gives labelled blocks of violins
#   with gaps, like facets (e.g. group.by = "condition", split.by =
#   "cell_type"), each marked by a bar and its name (R/group_bars.R), by
#   default opposite the group names (`split.by.side`). Seurat's dodged
#   violins are not drawn.
# * Orientation: `flip = FALSE` draws vertical violins (groups on x);
#   `flip = TRUE` horizontal ones (groups on y, the first at the top). The
#   meaning of `flip` is that of DotPlot_bioToolkit(); Seurat's `flip` only
#   applies to stacked violins. Ridges are horizontal only (ggridges).
# * `stack = TRUE`: one thin panel per feature, sharing the group axis (the
#   compact layout for many markers, as scanpy's stacked violin).
# * Sizes (F2): one violin is `cell.width` (vertical) or `cell.height`
#   (horizontal) inches wide; the value axis is `panel.height` /
#   `panel.width` inches long.
# * No p-values: tests on cells treat cells as independent replicates and
#   overstate significance.


#' @title Violin plot of genes or per-cell values
#'
#' @description Draws the distribution of one or more features in each group
#'              of cells as violins, with a box of the inter-quartile range
#'              and median and a sample of the cells as points. Features can
#'              be genes, metadata columns (QC metrics, module scores) or a
#'              mix, each in its own panel with its own axis. A drop-in
#'              alternative to \code{Seurat::VlnPlot()} with a fixed size:
#'              save it with \code{\link{save_plot}}; the values behind it are
#'              returned by \code{\link{export_source_data}}.
#'
#' @param object A Seurat object, or a data frame with one row per cell (e.g.
#'               a table of QC metrics; row names are the cell names)
#'
#' @param features Features to plot: genes, metadata columns, module scores
#'                 or embedding dimensions (from a data frame: its numeric
#'                 columns). Gene symbols are matched exactly (rat and mouse
#'                 "Cd3e", human "CD3E"); missing ones are left out with a
#'                 warning
#'
#' @param group.by Column to group cells by: one violin per group, filled
#'                 with the group's colour. If NULL (default), the active
#'                 identities of a Seurat object, or one group for a data frame
#'
#' @param split.by Column of a higher-level grouping, drawn like facets: one
#'                 labelled block of violins per level, e.g. `group.by =
#'                 "condition", split.by = "cell_type"` compares conditions
#'                 within each cell type. Default is NULL. (Seurat's
#'                 `split.by` dodges violins instead)
#'
#' @param split.by.cols Colours of the bars marking the `split.by` blocks:
#'                      NULL (default) for black bars, a qualitative palette
#'                      name, or colours named by level of `split.by`
#'
#' @param split.by.side Side of the `split.by` blocks: "opposite" (default;
#'                      opposite the group names: on top, or on the right with
#'                      `flip = TRUE`) or "labels" (beyond the group names)
#'
#' @param group.bar.width Thickness of the `split.by` bars in points, default
#'                        is 3.5
#'
#' @param cells Cells to use, as names or indices. If NULL (default), all
#'
#' @param assay Assay to read genes from. If NULL (default), the default
#'              assay; the assay and layer used are reported in a message
#'              when a feature is a gene, and written to the source data
#'
#' @param slot Layer to read genes from, default is "data"
#'
#' @param flip Whether to draw horizontal violins (groups on the y-axis, the
#'             first at the top), default is FALSE (vertical, groups on x)
#'
#' @param stack Whether to stack the features as thin panels sharing the
#'              group axis (rows; columns when flipped), default is FALSE
#'              (one panel each, in a grid of `ncol` columns)
#'
#' @param log Whether to draw values on a log10 axis: TRUE / FALSE for all
#'            features, or the names of the features to draw so (e.g.
#'            `c("nCount_RNA", "nFeature_RNA")`). Values of 0 or below cannot
#'            be drawn and are left out, with a message. Default is FALSE
#'
#' @param percent Whether values are fractions (0 to 1) shown as
#'                percentages on the axis: TRUE / FALSE, or the names of the
#'                features (e.g. "mito_fraction"). Default is FALSE
#'
#' @param y.max Top of the value axis: a number, a quantile of each feature's
#'              values such as "q99" (zooms past a few extreme cells), or
#'              NULL (default) for the full range. Only the view is cut:
#'              violins and boxes are computed from all cells and cut at the
#'              top of the axis, and points above it are not drawn. One value
#'              for all features, or a vector named by feature
#'
#' @param zeros How cells with a value of 0 are handled for genes: "include"
#'              (default, as Seurat), "exclude" (violins of the cells above
#'              0 only) or "annotate" (as "exclude", with the \% of cells
#'              above 0 printed for each group). Not applied to other
#'              features
#'
#' @param box Whether to draw the box (inter-quartile range and median),
#'            default is TRUE
#'
#' @param box.width Width of the box, as a fraction of a violin's slot,
#'                  default is 0.2
#'
#' @param pt.n Largest number of cells drawn as points per violin, picked at
#'             random (`seed`), default is 1500. 0 draws no points; Inf all
#'
#' @param pt.size Size of the points, default is 0.16
#'
#' @param alpha Opacity of the points, default is 0.4
#'
#' @param adjust Bandwidth multiplier of the densities, default is 1
#'
#' @param cols Colours of the groups: NULL (default) for the colours
#'             \code{\link{DimPlot_bioToolkit}} gives them ("scanpy-20" up to
#'             20 groups), a qualitative palette name, or colours named by
#'             group (or in level order)
#'
#' @param italic.genes Whether gene names are set in italics, default is TRUE
#'
#' @param ncol Number of columns of panels when `stack = FALSE`. If NULL
#'             (default), as Seurat
#'
#' @param cell.width,cell.height Width (vertical violins) or height
#'                               (horizontal violins) of one violin's slot in
#'                               inches, default is 0.35
#'
#' @param panel.width,panel.height Length of the value axis in inches: the
#'                                 panel height for vertical violins and the
#'                                 width for horizontal ones. If NULL
#'                                 (default), 1.6, or 0.5 when stacked
#'
#' @param raster Whether to rasterize the points: NULL (default) does so
#'               above 2000 points; or TRUE / FALSE
#'
#' @param raster.dpi Resolution of rasterized points, default is 1000
#'
#' @param seed Random seed of the sampled and jittered points, default is 42
#'
#' @param legend.position Position of the group legend, default is "none"
#'                        (the groups are named on the axis)
#'
#' @param bg Background colour, default is "transparent"
#'
#' @details
#' Cells without a value of `group.by` (or `split.by`) are left out, with a
#' message. Group colours are taken from all levels of `group.by` in the
#' whole object, so that they match an embedding of all cells.
#'
#' The source data (\code{\link{export_source_data}}) hold one row per cell
#' and feature: `cell`, `group`, `split`, `feature`, `value`, whether the cell
#' is `drawn` (`zeros`, `log`) and, for genes, the `assay` and `layer`.
#'
#' @return A ggplot object for one feature, otherwise a patchwork
#'
#' @examples
#' \dontrun{
#' # QC metrics per sample, from the metadata alone
#' VlnPlot_bioToolkit(obj[[]], c("nCount_RNA", "nFeature_RNA", "percent.mt"),
#'                    group.by = "sample", log = c("nCount_RNA", "nFeature_RNA"))
#'
#' # Genes per cell type, horizontal
#' VlnPlot_bioToolkit(obj, c("Cdh5", "Col1a1"), group.by = "cell_type",
#'                    assay = "SCT", flip = TRUE)
#'
#' # Conditions compared within each cell type; % above 0 for sparse genes
#' VlnPlot_bioToolkit(obj, "Nppa", group.by = "condition",
#'                    split.by = "cell_type", zeros = "annotate")
#'
#' # Many markers, stacked
#' p <- VlnPlot_bioToolkit(obj, c("Cdh5", "Col1a1", "Myh6", "Cd74"),
#'                         group.by = "cell_type", stack = TRUE, pt.n = 0)
#' save_plot(p, "violins.pdf")
#' }
#'
#' @seealso \code{\link{RidgePlot_bioToolkit}}, \code{\link{DotPlot_bioToolkit}},
#'          \code{\link{save_plot}}, \code{\link{export_source_data}}
#'
#' @import ggplot2
#' @export
#'
VlnPlot_bioToolkit <- function(
    object,
    features,
    group.by = NULL,
    split.by = NULL,
    split.by.cols = NULL,
    split.by.side = c("opposite", "labels"),
    group.bar.width = 3.5,
    cells = NULL,
    assay = NULL,
    slot = "data",
    flip = FALSE,
    stack = FALSE,
    log = FALSE,
    percent = FALSE,
    y.max = NULL,
    zeros = c("include", "exclude", "annotate"),
    box = TRUE,
    box.width = 0.2,
    pt.n = 1500,
    pt.size = 0.16,
    alpha = 0.4,
    adjust = 1,
    cols = NULL,
    italic.genes = TRUE,
    ncol = NULL,
    cell.width = 0.35,
    cell.height = 0.35,
    panel.width = NULL,
    panel.height = NULL,
    raster = NULL,
    raster.dpi = 1000,
    seed = 42,
    legend.position = "none",
    bg = "transparent"
) {
    zeros <- match.arg(zeros)
    split.by.side <- match.arg(split.by.side)
    check_bar_width(group.bar.width)
    if (!(is.numeric(pt.n) && length(pt.n) == 1 && pt.n >= 0)) {
        stop("`pt.n` must be a single number >= 0", call. = FALSE)
    }
    dd <- distribution_data(object, features, group.by, split.by, cells, assay, slot, zeros, log)
    features <- dd$features
    log <- feature_flags(log, features, "log")
    percent <- feature_flags(percent, features, "percent")
    group_cols <- group_colours(cols, dd$group_levels)$cols
    sz <- theme_sizes()
    # Default length of the value axis (inches)
    value_len <- if (stack) 0.5 else 1.6
    if (!flip && is.null(panel.height)) panel.height <- value_len
    if (flip && is.null(panel.width)) panel.width <- value_len

    # Points: at most pt.n cells per violin (group, and split), the same cells
    # for every feature
    pts_cells <- if (pt.n > 0) sample_per_unit(dd$cells, pt.n, seed) else character(0)
    n_points <- length(pts_cells) * length(features)
    layer_fun <- points_layer_fun(raster, n_points, 2000, raster.dpi)
    # Side of the split.by blocks (bars and names): opposite the group names
    # (top, or right when flipped) by default
    split_side <- if (!flip) {
        if (split.by.side == "opposite") "top" else "bottom"
    } else {
        if (split.by.side == "opposite") "right" else "left"
    }
    split_cols <- if (!is.null(dd$long$split)) {
        group_bar_colours(split.by.cols, levels(droplevels(as.factor(dd$long$split))), "split.by.cols")
    }

    plots <- lapply(seq_along(features), function(i) {
        f <- features[i]
        d <- dd$long[dd$long$feature == f, , drop = FALSE]
        first <- i == 1
        last <- i == length(features)
        # Top of the value axis (y.max), on the scale it is drawn on
        top <- value_limit(y.max, f, d$value[d$drawn])
        pts <- d[d$drawn & d$cell %in% pts_cells, , drop = FALSE]
        if (!is.na(top)) pts <- pts[pts$value <= top, , drop = FALSE]
        # Stacked panels show the split.by blocks once, on the panel at their
        # side (the first or last); the others have no block names (each
        # side's element is set, so each is blanked)
        show_split <- !stack || (if (split_side %in% c("top", "left")) first else last)
        no_strip_text <- theme(
            strip.text.x.top = element_blank(),
            strip.text.x.bottom = element_blank(),
            strip.text.y.left = element_blank(),
            strip.text.y.right = element_blank()
        )
        p <- violin_panel(
            d = d,
            pts = pts,
            group_cols = group_cols,
            flip = flip,
            box = box,
            box.width = box.width,
            pt.size = pt.size,
            alpha = alpha,
            adjust = adjust,
            layer_fun = layer_fun,
            seed = seed,
            annotate = zeros == "annotate" && dd$is_gene[[f]],
            top = if (is.na(top)) NA_real_ else if (log[[f]]) log10(top) else top,
            split_side = split_side
        )
        title <- feature_title(f, italic.genes && dd$is_gene[[f]], bold = !stack)
        p <- p +
            value_scale(if (flip) "x" else "y", log[[f]], percent[[f]], n.breaks = if (stack) 3 else 5,
                        far = if (zeros == "annotate" && dd$is_gene[[f]]) 0.2 else 0.05) +
            value_coord(flip, top) +
            group_axis_theme(flip, sz, group.bar.width) +
            theme(legend.position = legend.position)
        p <- p + distribution_size_theme(d, flip, cell.width, cell.height, panel.width, panel.height) +
            block_label_theme(d, flip, cell.width, sz, split_side, group.bar.width)
        if (!is.null(d$split) && show_split) {
            p <- p + split_bar_layer(d, flip, split_side, split_cols, group.bar.width, sz)
        }
        # Margins: none between stacked panels (they touch), the theme's
        # otherwise; more on the left when the first group label is long
        margins <- if (stack) c(1, 2, 1, 2) else rep(0.2 * 72 / 2.54, 4)
        if (!flip) margins <- label_overflow_margin(levels(droplevels(d$group))[1], cell.width, sz, margins)
        p <- p + theme(plot.margin = margin(margins[1], margins[2], margins[3], margins[4], unit = "pt"))
        value_title <- if (dd$is_gene[[f]] && !stack) "Expression" else NULL
        if (!stack) {
            p <- p + labs(title = title, x = if (flip) value_title else NULL, y = if (flip) NULL else value_title)
        } else if (!flip) {
            # Rows of features: the name to the left, set horizontally; group
            # names under the last row, block names above the first
            p <- p + labs(x = NULL, y = title) +
                theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = sz$base))
            if (!last) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
            if (!show_split) p <- p + no_strip_text
        } else {
            # Columns of features: the name above, set upright; group names
            # left of the first column, block names left of it too
            p <- p + labs(x = NULL, y = NULL, title = title) +
                theme(plot.title = element_text(angle = 90, hjust = 0, vjust = 0.5, size = sz$base,
                                                face = "plain"))
            if (!first) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
            if (!show_split) p <- p + no_strip_text
        }
        p <- p + theme_bg(bg)
        attr(p, "source_data") <- distribution_source(d, dd$assay)
        p
    })
    names(plots) <- features
    combine_feature_plots(plots, ncol, stack, flip, bg)
}


#' @title Ridge plot of genes or per-cell values
#'
#' @description Draws the distribution of one or more features in each group
#'              of cells as overlapping density ridges, the first group at the
#'              top. Features can be genes, metadata columns or a mix, each in
#'              its own panel with its own axis. A drop-in alternative to
#'              \code{Seurat::RidgePlot()} with a fixed size; it takes the same
#'              arguments as \code{\link{VlnPlot_bioToolkit}} where they apply.
#'              Needs the ggridges package.
#'
#' @inheritParams VlnPlot_bioToolkit
#'
#' @param group.by Column to group cells by: one ridge per group, filled with
#'                 the group's colour. If NULL (default), the active
#'                 identities of a Seurat object, or one group for a data frame
#'
#' @param y.max Right end of the value axis (Seurat's name; the value axis
#'              of a ridge plot is the x-axis): a number, a quantile of each
#'              feature's values such as "q99", or NULL (default) for the full
#'              range. One value for all features, or a vector named by
#'              feature
#'
#' @param zeros How cells with a value of 0 are handled for genes: "include"
#'              (default), "exclude" (ridges of the cells above 0 only) or
#'              "annotate" (as "exclude", with the \% of cells above 0
#'              printed at the right end of each ridge)
#'
#' @param overlap Height of the tallest ridge, in rows: above 1, ridges
#'                overlap the row above. Default is 1.8
#'
#' @param cell.height Distance between ridges (one row) in inches, default is
#'                    0.3
#'
#' @param panel.width Length of the value axis in inches, default is 2.5
#'
#' @details
#' Densities are computed with one bandwidth per feature (Silverman's rule on
#' all plotted values, times `adjust`), so that ridges of a feature are
#' smoothed alike. Each ridge is scaled to its own height, so the ridges
#' compare the shapes of the distributions, not the number of cells. The
#' source data are as in \code{\link{VlnPlot_bioToolkit}}.
#'
#' @return A ggplot object for one feature, otherwise a patchwork
#'
#' @examples
#' \dontrun{
#' RidgePlot_bioToolkit(obj[[]], c("nCount_RNA", "percent.mt"), group.by = "sample",
#'                      log = "nCount_RNA")
#' RidgePlot_bioToolkit(obj, "Col1a1", group.by = "cell_type", zeros = "annotate")
#' }
#'
#' @seealso \code{\link{VlnPlot_bioToolkit}}, \code{\link{save_plot}}
#'
#' @import ggplot2
#' @export
#'
RidgePlot_bioToolkit <- function(
    object,
    features,
    group.by = NULL,
    cells = NULL,
    assay = NULL,
    slot = "data",
    log = FALSE,
    percent = FALSE,
    y.max = NULL,
    zeros = c("include", "exclude", "annotate"),
    overlap = 1.8,
    adjust = 1,
    cols = NULL,
    italic.genes = TRUE,
    ncol = NULL,
    cell.height = 0.3,
    panel.width = 2.5,
    legend.position = "none",
    bg = "transparent"
) {
    check_installed("ggridges", reason = "to draw ridges")
    zeros <- match.arg(zeros)
    if (!(is.numeric(overlap) && length(overlap) == 1 && overlap > 0)) {
        stop("`overlap` must be a single positive number", call. = FALSE)
    }
    dd <- distribution_data(object, features, group.by, NULL, cells, assay, slot, zeros, log)
    features <- dd$features
    log <- feature_flags(log, features, "log")
    percent <- feature_flags(percent, features, "percent")
    group_cols <- group_colours(cols, dd$group_levels)$cols
    sz <- theme_sizes()

    plots <- lapply(features, function(f) {
        d <- dd$long[dd$long$feature == f, , drop = FALSE]
        annotate <- zeros == "annotate" && dd$is_gene[[f]]
        x_top <- value_limit(y.max, f, d$value[d$drawn])
        p <- ridge_panel(d, group_cols, log[[f]], adjust, overlap, annotate = annotate)
        rows <- attr(p, "ridge_rows")
        n <- length(rows$labels)
        # Rows are 1..n from the bottom; below the bottom ridge a fifth of a
        # row, above the top ridge its height plus a margin (at least half a
        # row, for the label)
        low <- 0.2
        high <- max(rows$top + 0.15, 0.5)
        p <- p +
            value_scale("x", log[[f]], percent[[f]], n.breaks = 5, far = if (annotate) 0.15 else 0.05) +
            scale_y_continuous(breaks = seq_len(n), labels = rows$labels, limits = c(1 - low, n + high),
                               expand = expansion(0)) +
            value_coord(TRUE, x_top) +
            theme_bioToolkit("axis", bg = bg) +
            panel_size_theme(panel.width, (n - 1 + low + high) * cell.height) +
            labs(title = feature_title(f, italic.genes && dd$is_gene[[f]], bold = TRUE), x = NULL, y = NULL) +
            theme(axis.line.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.text.y = element_text(vjust = 0),
                  legend.position = legend.position)
        attr(p, "source_data") <- distribution_source(d, dd$assay)
        p
    })
    names(plots) <- features
    combine_feature_plots(plots, ncol, stack = FALSE, flip = TRUE, bg)
}


#' Values of features per cell, with their groups, in long form
#'
#' @param object A Seurat object, or a data frame with one row per cell
#'
#' @param features,group.by,split.by,cells,assay,slot,zeros,log As in
#'        \code{\link{VlnPlot_bioToolkit}}
#'
#' @return A list with `long` (one row per cell and feature: `cell`,
#'         `group`, `split`, `feature`, `value`, `drawn`), `features` (those
#'         found), `is_gene`, `assay` (or NULL), `cells` (a data frame of the
#'         cells with their `unit`) and `group_levels` (a factor with all
#'         levels of `group.by`, for the colours)
#'
#' @keywords internal
distribution_data <- function(object,
                              features,
                              group.by,
                              split.by,
                              cells,
                              assay,
                              slot,
                              zeros,
                              log) {
    if (!(is.data.frame(object) || inherits(object, "Seurat"))) {
        stop("`object` must be a Seurat object or a data frame", call. = FALSE)
    }
    if (inherits(object, "Seurat")) check_installed("Seurat", reason = "to fetch values")
    cells <- plot_cells(object, cells)
    if (length(cells) == 0) stop("No cells to plot", call. = FALSE)
    all_cells <- plot_cells(object, NULL)
    group_all <- fetch_plot_group(object, group.by, all_cells)
    group <- group_all[match(cells, all_cells)]
    split <- NULL
    if (!is.null(split.by)) {
        if (length(split.by) != 1) stop("`split.by` must be a single column", call. = FALSE)
        split <- fetch_plot_group(object, split.by, cells)
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
    if (length(cells) == 0) stop("No cells to plot", call. = FALSE)

    fv <- fetch_plot_values(object, features, cells, assay, slot)
    features <- colnames(fv$data)
    log_f <- feature_flags(log, features, "log")
    group <- droplevels(group)
    if (!is.null(split)) split <- droplevels(split)

    long <- do.call(rbind, lapply(features, function(f) {
        v <- fv$data[[f]]
        drawn <- !is.na(v)
        # Sparse genes: cells at 0 are left out of the violins / ridges
        if (zeros != "include" && fv$is_gene[[f]]) drawn <- drawn & v > 0
        if (log_f[[f]]) {
            n_bad <- sum(drawn & v <= 0)
            if (n_bad > 0) message(n_bad, " cell(s) with '", f, "' <= 0 cannot be drawn on a log axis")
            drawn <- drawn & v > 0
        }
        data.frame(cell = cells, group = group, split = if (is.null(split)) NA else split,
                   feature = f, value = v, drawn = drawn, stringsAsFactors = FALSE)
    }))
    if (is.null(split)) long$split <- NULL
    unit <- if (is.null(split)) group else interaction(split, group, drop = TRUE, lex.order = TRUE)
    list(
        long = long,
        features = features,
        is_gene = fv$is_gene,
        assay = if (!is.null(fv$assay)) list(assay = fv$assay$assay, layer = slot),
        cells = data.frame(cell = cells, unit = unit, stringsAsFactors = FALSE),
        group_levels = group_all
    )
}


#' Resolve a per-feature switch given as TRUE / FALSE or feature names
#'
#' @param x TRUE, FALSE, or names of features
#'
#' @param features All features
#'
#' @param arg Name of the argument, for messages
#'
#' @return A logical vector named by feature
#'
#' @keywords internal
feature_flags <- function(x,
                          features,
                          arg) {
    if (is.logical(x) && length(x) == 1 && !is.na(x)) {
        return(stats::setNames(rep(x, length(features)), features))
    }
    if (!is.character(x)) stop("`", arg, "` must be TRUE, FALSE or feature names", call. = FALSE)
    unknown <- setdiff(x, features)
    if (length(unknown)) {
        warning("`", arg, "` names feature(s) not plotted: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    stats::setNames(features %in% x, features)
}


#' Pick at most `n` cells per unit at random
#'
#' @param cells A data frame with `cell` and `unit`
#'
#' @param n Largest number of cells per unit (Inf for all)
#'
#' @param seed Random seed
#'
#' @return Cell names
#'
#' @keywords internal
sample_per_unit <- function(cells,
                            n,
                            seed) {
    by_unit <- split(cells$cell, cells$unit, drop = TRUE)
    with_seed(seed, unlist(lapply(by_unit, function(x) {
        if (length(x) <= n) x else x[sample.int(length(x), n)]
    }), use.names = FALSE))
}


#' Title of a feature panel, in italics for genes
#'
#' @param f Feature name
#'
#' @param italic Whether to set it in italics
#'
#' @param bold Whether to set it in bold
#'
#' @return A string or a plotmath expression
#'
#' @keywords internal
feature_title <- function(f,
                          italic,
                          bold) {
    if (!italic) return(f)
    if (bold) bquote(bolditalic(.(f))) else bquote(italic(.(f)))
}


#' Violins, points and boxes of one feature
#'
#' @param d Long table of one feature (\code{distribution_data()})
#'
#' @param pts Rows of `d` drawn as points
#'
#' @param group_cols Colours named by group
#'
#' @param flip Whether violins are horizontal
#'
#' @param box,box.width,pt.size,alpha,adjust,seed As in
#'        \code{\link{VlnPlot_bioToolkit}}
#'
#' @param layer_fun Function rasterizing the point layer, or not
#'
#' @param annotate Whether to print the \% of cells above 0 per group
#'
#' @param top Top of the value axis on the scale it is drawn on (log10 for a
#'            log axis), where violins and boxes are cut; NA for none
#'
#' @param split_side Side of the split.by block names: "top" / "bottom"
#'                   (vertical violins) or "left" / "right" (horizontal)
#'
#' @return A ggplot object (no scales for the value axis yet)
#'
#' @keywords internal
violin_panel <- function(d,
                         pts,
                         group_cols,
                         flip,
                         box,
                         box.width,
                         pt.size,
                         alpha,
                         adjust,
                         layer_fun,
                         seed,
                         annotate,
                         top = NA_real_,
                         split_side = if (flip) "right" else "top") {
    lw <- pt_to_linewidth(theme_sizes()$line)
    # Horizontal violins list the first group at the top
    group_lv <- levels(droplevels(d$group))
    lv <- if (flip) rev(group_lv) else group_lv
    fix <- function(x) {
        x$group <- factor(as.character(x$group), levels = lv)
        x
    }
    dv <- fix(d[d$drawn, , drop = FALSE])
    pts <- fix(pts)
    map <- if (flip) aes(x = .data$value, y = .data$group) else aes(x = .data$group, y = .data$value)
    # Every group keeps its slot, even one with no cell drawn (e.g. no cell
    # above 0), at a value that is drawn (so that it is valid on a log axis)
    blank <- fix(unique(d[c("group", if (!is.null(d$split)) "split")]))
    blank$value <- if (nrow(dv) > 0) dv$value[1] else 1
    # A density needs two cells: groups with fewer get their box and points
    # only
    n_drawn <- stats::ave(seq_len(nrow(dv)), dv$group, if (is.null(dv$split)) 1 else dv$split, FUN = length)
    # Violins (and boxes) are cut at `top` by their geoms, the panel being
    # unclipped
    violins <- geom_violin(
        aes(fill = .data$group),
        data = dv[n_drawn >= 2, , drop = FALSE],
        scale = "width",
        width = 0.88,
        colour = NA,
        alpha = 0.9,
        trim = TRUE,
        adjust = adjust,
        orientation = if (flip) "y" else "x"
    )
    violins$geom <- GeomViolinTop
    violins$geom_params$top <- top
    p <- ggplot(dv, map) +
        geom_blank(data = blank) +
        violins
    if (nrow(pts) > 0) {
        jitter <- if (flip) {
            position_jitter(width = 0, height = 0.13, seed = seed)
        } else {
            position_jitter(width = 0.13, height = 0, seed = seed)
        }
        p <- p + layer_fun(geom_point(data = pts, position = jitter, shape = 16, stroke = 0,
                                      size = pt.size, alpha = alpha, colour = "grey15"))
    }
    if (box) {
        # coef = 0: no whiskers, the box only (quartiles and median)
        boxes <- geom_boxplot(
            width = box.width,
            outlier.shape = NA,
            fill = "white",
            colour = "grey15",
            linewidth = lw,
            coef = 0,
            orientation = if (flip) "y" else "x"
        )
        boxes$geom <- GeomBoxplotTop
        boxes$geom_params$top <- top
        p <- p + boxes
    }
    if (annotate) {
        keys <- c("group", if (!is.null(d$split)) "split")
        pct <- stats::aggregate(list(pct = d$value > 0 & !is.na(d$value)), d[keys], mean)
        pct <- fix(pct)
        pct$label <- pct_label(pct$pct)
        text_size <- pt_to_size(theme_sizes()$small * 0.85)
        p <- p + if (flip) {
            geom_text(data = pct, aes(x = Inf, y = .data$group, label = .data$label), inherit.aes = FALSE,
                      hjust = 1.1, size = text_size, colour = "grey30")
        } else {
            geom_text(data = pct, aes(x = .data$group, y = Inf, label = .data$label), inherit.aes = FALSE,
                      vjust = 1.3, size = text_size, colour = "grey30")
        }
    }
    p <- p + scale_fill_manual(values = group_cols, drop = TRUE, guide = guide_legend(reverse = flip))
    if (!is.null(d$split)) {
        # Block names on top / right by default; at the bottom / left on the
        # side of the group names
        p <- p + if (flip) {
            facet_grid(
                rows = vars(.data$split),
                scales = "free_y",
                space = "free_y",
                switch = if (split_side == "left") "y"
            )
        } else {
            facet_grid(
                cols = vars(.data$split),
                scales = "free_x",
                space = "free_x",
                switch = if (split_side == "bottom") "x"
            )
        }
    }
    p <- p + if (flip) {
        scale_y_discrete(expand = grid_expansion())
    } else {
        scale_x_discrete(expand = grid_expansion())
    }
    p
}


#' Ridges of one feature
#'
#' @param d Long table of one feature (\code{distribution_data()})
#'
#' @param group_cols Colours named by group
#'
#' @param log Whether the value axis is log10 (the bandwidth is then chosen
#'            on that scale)
#'
#' @param adjust Bandwidth multiplier
#'
#' @param overlap Height of the tallest ridge, in rows
#'
#' @param annotate Whether to print the \% of cells above 0 per group
#'
#' @return A ggplot object (no scales for the value axis yet)
#'
#' @keywords internal
ridge_panel <- function(d,
                        group_cols,
                        log,
                        adjust,
                        overlap,
                        annotate) {
    # Rows are numbered from the bottom: the first group at the top. Rows
    # are numeric (not a discrete scale), so that a group with no ridge (e.g.
    # no cell above 0) keeps its row and the others are not shifted
    lv <- rev(levels(droplevels(d$group)))
    dv <- d[d$drawn, , drop = FALSE]
    # A density needs two cells
    n_drawn <- table(factor(dv$group, levels = lv))
    dv <- dv[dv$group %in% names(n_drawn)[n_drawn >= 2], , drop = FALSE]
    dv$pos <- match(as.character(dv$group), lv)
    dv$group <- factor(as.character(dv$group), levels = lv)
    v <- if (log) log10(dv$value) else dv$value
    # One bandwidth for the whole feature, as ggridges' "joint" bandwidth
    # (which it otherwise announces in a message), times `adjust`
    bw <- if (length(unique(v)) > 1) stats::bw.nrd0(v) * adjust else 1
    p <- ggplot(dv, aes(x = .data$value, y = .data$pos, group = .data$group, fill = .data$group)) +
        ggridges::geom_density_ridges(scale = overlap, rel_min_height = 0.004, bandwidth = bw,
                                      colour = "white", linewidth = pt_to_linewidth(0.6), alpha = 0.92) +
        scale_fill_manual(values = group_cols, drop = TRUE, guide = guide_legend(reverse = TRUE))
    # Height of the top ridge, in rows: ggridges scales the tallest ridge of
    # the panel to `overlap` rows, so the top one reaches overlap times its
    # peak density relative to the tallest
    peaks <- vapply(split(v, droplevels(dv$group)), function(x) {
        if (length(unique(x)) < 2) return(0)
        max(stats::density(x, bw = bw, from = min(v), to = max(v))$y)
    }, numeric(1))
    top <- if (length(peaks) && lv[length(lv)] %in% names(peaks) && max(peaks) > 0) {
        overlap * peaks[[lv[length(lv)]]] / max(peaks)
    } else {
        0
    }
    attr(p, "ridge_rows") <- list(labels = lv, top = top)
    if (annotate) {
        pct <- stats::aggregate(list(pct = d$value > 0 & !is.na(d$value)), d["group"], mean)
        pct$pos <- match(as.character(pct$group), lv)
        pct$label <- pct_label(pct$pct)
        p <- p + geom_text(data = pct, aes(x = Inf, y = .data$pos, label = .data$label),
                           inherit.aes = FALSE, hjust = 1, vjust = -0.4,
                           size = pt_to_size(theme_sizes()$small * 0.85), colour = "grey30")
    }
    p
}


#' Continuous scale of the value axis
#'
#' @param axis "x" or "y"
#'
#' @param log Whether log10
#'
#' @param percent Whether values are fractions labelled as percentages
#'
#' @param n.breaks Desired number of breaks
#'
#' @param far Expansion of the high end, as a fraction of the range (room
#'            for labels there)
#'
#' @return A ggplot2 scale
#'
#' @keywords internal
value_scale <- function(axis,
                        log,
                        percent,
                        n.breaks,
                        far = 0.05) {
    labels <- if (percent) {
        scales::label_percent(drop0trailing = TRUE)
    } else {
        # Large counts as 10K, 1.5M
        scales::label_number(scale_cut = scales::cut_short_scale(), drop0trailing = TRUE)
    }
    fun <- if (axis == "x") {
        if (log) scale_x_log10 else scale_x_continuous
    } else {
        if (log) scale_y_log10 else scale_y_continuous
    }
    expand <- expansion(mult = c(0.05, far))
    if (log) fun(labels = labels, expand = expand) else fun(labels = labels, n.breaks = n.breaks, expand = expand)
}


#' Top of the value axis from `y.max`
#'
#' @param y.max NULL, a number or quantile ("q99"), or a vector named by
#'              feature
#'
#' @param f Feature name
#'
#' @param v Values of the feature that are drawn
#'
#' @return A number, or NA for the full range
#'
#' @keywords internal
value_limit <- function(y.max,
                        f,
                        v) {
    if (is.null(y.max)) return(NA_real_)
    if (!is.null(names(y.max))) {
        if (!(f %in% names(y.max))) return(NA_real_)
        y.max <- y.max[[f]]
    } else if (length(y.max) != 1) {
        stop("`y.max` must be one value, or values named by feature", call. = FALSE)
    }
    cutoff_value(y.max, v[is.finite(v)], function(x) NA_real_)
}


#' Coordinates of a violin or ridge panel, zoomed to `y.max`
#'
#' @param flip Whether the value axis is x
#'
#' @param top Top of the value axis, or NA
#'
#' @return A coord
#'
#' @keywords internal
value_coord <- function(flip,
                        top,
                        clip = "off") {
    lim <- if (is.na(top)) NULL else c(NA, top)
    if (flip) coord_cartesian(xlim = lim, clip = clip) else coord_cartesian(ylim = lim, clip = clip)
}


#' Cut a density (violin) at the top of the value axis
#'
#' @param d Rows of one violin after the stat (`y` sorted or not, with
#'          `violinwidth`), in the orientation of `flip_data()`
#'
#' @param top Top of the axis, on the scale drawn
#'
#' @return The rows up to `top`, closed by a row at `top` with the width
#'         interpolated there
#'
#' @keywords internal
cut_density <- function(d,
                        top) {
    d <- d[order(d$y), , drop = FALSE]
    keep <- d$y <= top
    if (all(keep)) return(d)
    if (!any(keep)) return(d[0, , drop = FALSE])
    edge <- d[max(which(keep)), , drop = FALSE]
    edge$violinwidth <- stats::approx(d$y, d$violinwidth, xout = top, ties = mean)$y
    edge$y <- top
    rbind(d[keep, , drop = FALSE], edge)
}


# Violins and boxes cut at the top of the value axis (`top`, on the scale
# drawn; NA for none). The panels of a violin plot are not clipped (the
# split.by bars sit outside them), so without this a violin taller than a
# lowered `y.max` would run past the axis
GeomViolinTop <- ggproto("GeomViolinTop", GeomViolin,
    extra_params = c(GeomViolin$extra_params, "top"),
    setup_data = function(data, params) {
        data <- GeomViolin$setup_data(data, params)
        top <- params$top
        if (is.null(top) || is.na(top) || nrow(data) == 0) return(data)
        flipped <- isTRUE(data$flipped_aes[1])
        data <- flip_data(data, flipped)
        # One violin per group and panel (group ids repeat across facets)
        parts <- split(data, list(data$PANEL, data$group), drop = TRUE)
        data <- do.call(rbind, lapply(parts, cut_density, top = top))
        rownames(data) <- NULL
        flip_data(data, flipped)
    }
)

GeomBoxplotTop <- ggproto("GeomBoxplotTop", GeomBoxplot,
    extra_params = c(GeomBoxplot$extra_params, "top"),
    setup_data = function(data, params) {
        data <- GeomBoxplot$setup_data(data, params)
        top <- params$top
        if (is.null(top) || is.na(top) || nrow(data) == 0) return(data)
        flipped <- isTRUE(data$flipped_aes[1])
        data <- flip_data(data, flipped)
        # A box (or part of it) above the top is cut there; one wholly above
        # it is not drawn
        data <- data[data$lower <= top, , drop = FALSE]
        for (v in intersect(c("ymin", "lower", "middle", "upper", "ymax", "ymax_final"), names(data))) {
            data[[v]] <- pmin(data[[v]], top)
        }
        flip_data(data, flipped)
    }
)


#' Theme of the group axis of a violin plot
#'
#' @param flip Whether groups are on the y-axis
#'
#' @param sz Sizes from \code{theme_sizes()}
#'
#' @return A theme
#'
#' @keywords internal
group_axis_theme <- function(flip,
                             sz,
                             bar.width = 3.5) {
    block_text <- function(side, ...) {
        element_text(size = sz$small, face = "bold", margin = group_bar_margin(side, bar.width), ...)
    }
    th <- theme_bioToolkit("axis") +
        theme(
            # Block names beyond their bars (R/group_bars.R); on the side of
            # the group names, outside them
            strip.placement = "outside",
            strip.switch.pad.grid = unit(0, "pt"),
            strip.clip = "off",
            strip.text.x.top = block_text("top"),
            strip.text.x.bottom = block_text("bottom"),
            strip.text.y.left = block_text("left", angle = 0, hjust = 1),
            strip.text.y.right = block_text("right", angle = 0, hjust = 0),
            panel.spacing = unit(0.08, "in"),
            # Known label margins, so that a bar beyond the labels can be
            # placed (split_bar_layer())
            axis.text.x = element_text(margin = margin(t = 2)),
            axis.text.y = element_text(margin = margin(r = 2))
        )
    if (!flip) th <- th + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    th
}


#' Bars marking the split.by blocks of a violin plot
#'
#' @param d Long table of one feature
#'
#' @param flip Whether violins are horizontal
#'
#' @param side Side of the bars
#'
#' @param cols Bar colours named by block
#'
#' @param width Thickness of the bars in points
#'
#' @param sz Sizes from \code{theme_sizes()}
#'
#' @return A ggplot2 layer
#'
#' @keywords internal
split_bar_layer <- function(d,
                            flip,
                            side,
                            cols,
                            width,
                            sz) {
    d$split <- droplevels(as.factor(d$split))
    blocks <- levels(d$split)
    group_lv <- levels(droplevels(d$group))
    # Horizontal violins list the first group at the top (violin_panel())
    lv <- if (flip) rev(group_lv) else group_lv
    ends <- lapply(blocks, function(b) {
        g <- intersect(lv, as.character(d$group[d$split == b]))
        c(g[1], g[length(g)])
    })
    first <- factor(vapply(ends, `[`, "", 1), levels = lv)
    last <- factor(vapply(ends, `[`, "", 2), levels = lv)
    bars <- data.frame(split = factor(blocks, levels = blocks), fill = unname(cols[blocks]))
    if (flip) {
        bars$y <- first
        bars$yend <- last
    } else {
        bars$x <- first
        bars$xend <- last
    }
    gap <- grid::unit(GROUP_BAR_GAP, "pt")
    # On the side of the group names, beyond the ticks and the names
    offset <- if (side %in% c("top", "right")) {
        gap
    } else {
        labels <- unique(as.character(d$group))
        grid::unit(sz$small / 3 + 2, "pt") +
            label_extent(labels, sz$small, angle = if (flip) 0 else 45,
                         across = if (flip) "width" else "height") +
            gap
    }
    group_bar_layer(bars, side = side, width = width, offset = offset)
}


#' Label of the share of cells above 0
#'
#' @param x Fractions
#'
#' @return Labels such as "26\%", "<1\%" (above 0 but below 0.5\%) or "0\%"
#'
#' @keywords internal
pct_label <- function(x) {
    out <- paste0(format(round(x * 100), trim = TRUE), "%")
    out[x > 0 & x < 0.005] <- "<1%"
    out
}


#' Turn block (split) labels upright when they are wider than their block
#'
#' Widths are estimated from the number of characters (bold text, about 0.6
#' em per character), as in the dot grid.
#'
#' @param d Long table of one feature
#'
#' @param flip Whether blocks are rows (labels to the left, never turned)
#'
#' @param cell.width Width of one violin's slot
#'
#' @param sz Sizes from \code{theme_sizes()}
#'
#' @param side Side of the block names ("top" or "bottom" when not flipped)
#'
#' @param bar.width Thickness of the block bars in points
#'
#' @return A theme
#'
#' @keywords internal
block_label_theme <- function(d,
                              flip,
                              cell.width,
                              sz,
                              side = "top",
                              bar.width = 3.5) {
    if (flip || is.null(d$split)) return(theme())
    blocks <- split(as.character(d$group), droplevels(as.factor(d$split)))
    widths <- vapply(blocks, function(x) length(unique(x)), integer(1)) * cell.width
    label_w <- nchar(names(blocks)) * 0.6 * sz$small / 72
    if (!any(label_w > widths + 0.08)) return(theme())
    upright <- element_text(
        size = sz$small,
        face = "bold",
        angle = 90,
        hjust = if (side == "top") 0 else 1,
        margin = group_bar_margin(side, bar.width)
    )
    if (side == "top") theme(strip.text.x.top = upright) else theme(strip.text.x.bottom = upright)
}


#' Left margin for a first group label that runs past the plot's left edge
#'
#' Group labels of vertical violins are turned 45 degrees and end under their
#' violin, so a long first label reaches left of the panel, past the axis
#' and out of the saved figure. The margin is estimated from the number of
#' characters (about 0.55 em per character) less the room left of the
#' violin (half a slot and about 0.3 in for the axis).
#'
#' @param first First group label
#'
#' @param cell.width Width of one violin's slot
#'
#' @param sz Sizes from \code{theme_sizes()}
#'
#' @param margins Plot margins (top, right, bottom, left) in points
#'
#' @return The margins, with the left one widened if needed
#'
#' @keywords internal
label_overflow_margin <- function(first,
                                  cell.width,
                                  sz,
                                  margins) {
    reach <- nchar(first) * 0.55 * sz$small / 72 * cos(pi / 4)
    extra <- reach - cell.width / 2 - 0.3
    if (extra > 0) margins[4] <- margins[4] + extra * 72
    margins
}


#' Fixed panel sizes of a violin plot
#'
#' @param d Long table of one feature
#'
#' @param flip Whether violins are horizontal
#'
#' @param cell.width,cell.height,panel.width,panel.height As in
#'        \code{\link{VlnPlot_bioToolkit}}
#'
#' @return A theme
#'
#' @keywords internal
distribution_size_theme <- function(d,
                                    flip,
                                    cell.width,
                                    cell.height,
                                    panel.width,
                                    panel.height) {
    # Number of groups in each block (facet), as drawn with free scales
    n <- if (is.null(d$split)) {
        length(unique(d$group))
    } else {
        vapply(split(as.character(d$group), droplevels(as.factor(d$split))),
               function(x) length(unique(x)), integer(1))
    }
    if (!flip) {
        panel_size_theme(grid_panel_size(n, cell.width), panel.height)
    } else {
        panel_size_theme(panel.width, grid_panel_size(n, cell.height))
    }
}


#' Background of a plot
#'
#' @param bg Background colour
#'
#' @return A theme
#'
#' @keywords internal
theme_bg <- function(bg) {
    transparent <- is.na(bg) || identical(bg, "transparent")
    theme(plot.background = element_rect(fill = if (transparent) "transparent" else bg, colour = NA))
}


#' Source data of one feature panel
#'
#' @param d Long table of one feature
#'
#' @param assay NULL, or a list with `assay` and `layer`
#'
#' @return A data frame
#'
#' @keywords internal
distribution_source <- function(d,
                                assay) {
    if (!is.null(assay)) {
        d$assay <- assay$assay
        d$layer <- assay$layer
    }
    rownames(d) <- NULL
    d
}


#' Combine one plot per feature
#'
#' @param plots Named list of ggplot objects
#'
#' @param ncol Columns of the grid when not stacked, or NULL
#'
#' @param stack Whether to stack them as rows (columns when `flip`)
#'
#' @param flip Whether groups are on the y-axis
#'
#' @param bg Background colour
#'
#' @return The plot itself for one feature, otherwise a patchwork
#'
#' @keywords internal
combine_feature_plots <- function(plots,
                                  ncol,
                                  stack,
                                  flip,
                                  bg) {
    if (length(plots) == 1) return(plots[[1]])
    layout <- if (stack) {
        if (flip) patchwork::plot_layout(nrow = 1) else patchwork::plot_layout(ncol = 1)
    } else {
        patchwork::plot_layout(ncol = if (is.null(ncol)) seurat_feature_ncol(length(plots)) else ncol)
    }
    patchwork::wrap_plots(plots) + layout + patchwork::plot_annotation(theme = theme_bg(bg))
}
