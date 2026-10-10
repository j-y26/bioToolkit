# Functions to plot better reduction plots for Seurat objects

#' @title Publication-quality dimensionality reduction plot
#'
#' @description Plots cells on a 2D reduction (UMAP, tSNE, PCA, ...) coloured by
#'              a grouping variable, optionally split into one panel per level
#'              of a second variable. Every panel (the scatter itself) has a
#'              fixed absolute size set by `panel.width` and `panel.height`, so
#'              panels are identical across figures regardless of how many there
#'              are and of the legend's size or position. Split plots share a
#'              single legend covering all groups, even when a panel contains
#'              only some of them. Use \code{\link{save_dimplot}} to save it: the
#'              output file is sized to fit the panels, strips, legend and
#'              margins exactly.
#'
#' @param object A Seurat object
#'
#' @param group.by Metadata column to colour cells by. If NULL, the active
#'                 identities are used
#'
#' @param split.by Metadata column to split the plot by, one panel per level.
#'                 If NULL (default), a single panel is drawn. Levels with no
#'                 cells are dropped
#'
#' @param reduction Name of the reduction to plot, default is "umap"
#'
#' @param dims The two dimensions of the reduction to plot, default is 1 and 2
#'
#' @param cells Cells to plot, as names or indices (e.g. from
#'              \code{Seurat::WhichCells()}). Other cells are left out of the
#'              plot entirely. If NULL (default), all cells are plotted
#'
#' @param cells.highlight Cells to colour, as names or indices. Other plotted
#'                        cells are drawn in `na.value` and left out of the
#'                        legend. If NULL (default), all plotted cells are
#'                        coloured
#'
#' @param full.extent Whether the axes span all cells of the object rather than
#'                    only the plotted ones, so that a plot of a subset of cells
#'                    lines up with a plot of the whole object. Default is FALSE
#'
#' @param cols Colours for the groups, either named by group or in the order of
#'             the group levels. If NULL, the default ggplot2 hue palette is used
#'
#' @param na.value Colour of cells with no group and of cells not in
#'                 `cells.highlight`, default is "grey88"
#'
#' @param pt.size Point size, default is 0.15
#'
#' @param alpha Point transparency, default is 1
#'
#' @param shuffle Whether to plot cells in random order so that no group is
#'                systematically drawn on top, default is TRUE
#'
#' @param outline Whether to draw an outline around the cell clouds, default is
#'                FALSE
#'
#' @param outline.color Colour of the outline: "black" (default), "darken" for a
#'                      darkened version of each group's colour, or any colour
#'
#' @param outline.size Thickness of the outline (in point size units), default
#'                     is 0.6
#'
#' @param outline.darken Amount to darken group colours by when
#'                       `outline.color = "darken"`, between 0 and 1, default
#'                       is 0.4
#'
#' @param label Whether to label groups at their median position, default is
#'              FALSE. When split, groups are labelled within each panel
#'
#' @param label.size Size of the group labels, default is 3.5
#'
#' @param legend.position Position of the legend, one of "right" (default),
#'                        "left", "top", "bottom", "inside" or "none"
#'
#' @param legend.ncol Number of legend columns. If NULL (default), ggplot2
#'                    chooses based on the number of groups and the position
#'
#' @param legend.pt.size Size of the legend keys, default is 3.5
#'
#' @param legend.text.size Font size of the legend text, default is 11
#'
#' @param ncol Number of panel columns when `split.by` is given. If NULL
#'             (default), all panels are laid out in a single row
#'
#' @param by.row Whether panels fill the grid row by row (default, TRUE) or
#'               column by column (FALSE)
#'
#' @param panel.width Width of a single panel in inches, default is 4
#'
#' @param panel.height Height of a single panel in inches, default is
#'                     `panel.width`
#'
#' @param panel.spacing Gap between panels in inches, default is 0.1
#'
#' @param strip.text.size Font size of the panel (split) labels, default is 12
#'
#' @param strip.text.face Font face of the panel labels, default is "bold"
#'
#' @param raster Whether to rasterize the points while keeping text, legend
#'               and arrows as vectors: NULL (default) rasterizes when there
#'               are more than `raster.threshold` cells, or TRUE / FALSE
#'
#' @param raster.threshold Number of cells above which `raster = NULL`
#'                         rasterizes, default is 10000
#'
#' @param raster.dpi Resolution of the rasterized points, default is 1000. Match
#'                   the `dpi` used when saving raster (e.g. png) output
#'
#' @param bg Background colour, default is "white". Use "transparent" for a
#'           transparent background
#'
#' @param seed Random seed used for shuffling and label placement, default is
#'             42. The global random number generator state is not changed
#'
#' @param ... Additional arguments passed to \code{\link{add_dimension_arrows}}.
#'            `rasterize` is accepted as a deprecated alias of `raster`
#'
#' @return a ggplot object
#'
#' @import ggplot2
#' @export
#'
DimPlot_bioToolkit <- function(
    object,
    group.by = NULL,
    split.by = NULL,
    reduction = "umap",
    dims = c(1, 2),
    cells = NULL,
    cells.highlight = NULL,
    full.extent = FALSE,
    cols = NULL,
    na.value = "grey88",
    pt.size = 0.15,
    alpha = 1,
    shuffle = TRUE,
    outline = FALSE,
    outline.color = "black",
    outline.size = 0.6,
    outline.darken = 0.4,
    label = FALSE,
    label.size = 3.5,
    legend.position = c("right", "left", "top", "bottom", "inside", "none"),
    legend.ncol = NULL,
    legend.pt.size = 3.5,
    legend.text.size = 11,
    ncol = NULL,
    by.row = TRUE,
    panel.width = 4,
    panel.height = panel.width,
    panel.spacing = 0.1,
    strip.text.size = 12,
    strip.text.face = "bold",
    raster = NULL,
    raster.threshold = 1e4,
    raster.dpi = 1000,
    bg = "white",
    seed = 42,
    ...
) {

    # Preparation
    check_fixed_panel_support()
    legend.position <- match.arg(legend.position)
    arrow_args <- list(...)
    raster <- resolve_raster_arg(raster, arrow_args)
    arrow_args$rasterize <- NULL
    dims <- resolve_dims(dims, reduction)

    # Embedding of the plotted cells, with axis limits and split panels
    rd <- reduction_data(object, reduction, dims$dims, cells, split.by, full.extent)
    df <- rd$df

    # Fetch groups for the plotted cells so rows are guaranteed to align
    group <- if (is.null(group.by)) {
        Seurat::Idents(object)[rownames(df)]
    } else {
        fetch_meta(object, group.by, rownames(df))
    }
    df$group <- if (is.factor(group)) group else factor(group)
    df$highlight <- highlight_flags(object, rownames(df), cells.highlight)

    # Resolve colours up front so outlines can reuse them. Levels are kept for
    # all plotted cells, not per panel or highlight, so that every panel of a
    # split plot is coloured the same way and colours do not shift with the
    # highlighted cells
    if (is.null(cols)) {
        df$group <- droplevels(df$group)
        cols <- stats::setNames(scales::hue_pal()(nlevels(df$group)), levels(df$group))
    } else if (is.null(names(cols))) {
        if (length(cols) < nlevels(df$group)) {
            stop("`cols` has ", length(cols), " colours but there are ",
                 nlevels(df$group), " groups", call. = FALSE)
        }
        cols <- stats::setNames(cols[seq_len(nlevels(df$group))], levels(df$group))
    } else if (length(missing_cols <- setdiff(levels(df$group), names(cols)))) {
        stop("No colour given for group(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    # The legend lists the groups present among the coloured cells only
    legend_breaks <- levels(droplevels(df$group[df$highlight]))

    if (shuffle) df <- with_seed(seed, df[sample.int(nrow(df)), ])

    layout <- panel_layout(rd$n_panels, ncol, by.row, default_ncol = rd$n_panels)
    layer_fun <- points_layer_fun(raster, nrow(df), raster.threshold, raster.dpi)
    halo.color <- halo_colour(bg)

    p <- ggplot(df, aes(x = x, y = y))

    if (outline) {
        if (identical(outline.color, "darken")) {
            check_installed("colorspace", reason = "to darken outline colours")
            dark <- colorspace::darken(c(cols[levels(df$group)], na.value), amount = outline.darken)
            outline.color <- unname(dark[as.integer(df$group)])
            outline.color[is.na(df$group) | !df$highlight] <- dark[length(dark)]
        }
        p <- p + outline_layers(df, pt.size, outline.size, outline.color, halo.color, layer_fun)
    }

    # Cells outside the highlight are drawn first, in a constant colour that is
    # not mapped and therefore stays out of the legend
    if (any(!df$highlight)) {
        p <- p + layer_fun(geom_point(
            data = df[!df$highlight, ], colour = na.value,
            size = pt.size, alpha = alpha, stroke = 0
        ))
    }

    p <- p +
        layer_fun(geom_point(
            aes(colour = group), data = df[df$highlight, ],
            size = pt.size, alpha = alpha, stroke = 0
        )) +
        # drop = FALSE keeps every group in the legend even when a panel, or the
        # whole plot, does not contain all of them
        scale_colour_manual(values = cols, breaks = legend_breaks, na.value = na.value, drop = FALSE) +
        guides(colour = guide_legend(
            ncol = legend.ncol,
            override.aes = list(size = legend.pt.size, alpha = 1)
        )) +
        extent_scales(rd) +
        fixed_panel_theme(
            panel.width = panel.width, panel.height = panel.height,
            ncol = layout$ncol, nrow = layout$nrow,
            panel.spacing = panel.spacing, legend.position = legend.position,
            legend.text.size = legend.text.size, strip.text.size = strip.text.size,
            strip.text.face = strip.text.face, bg = bg
        ) +
        theme(
            legend.title = element_blank(),
            legend.key.size = unit(0.45, "cm"),
            legend.key.spacing.y = unit(0.05, "cm")
        )

    if (rd$split_plot) {
        p <- p + facet_wrap(~ panel, ncol = layout$ncol, nrow = layout$nrow, dir = layout$dir)
    }

    if (label) {
        p <- p + group_label_layer(
            df[df$highlight, ], by_panel = rd$split_plot, label.size = label.size,
            label.color = "black", halo.color = halo.color, seed = seed
        )
    }

    # Arrows in place of axes, drawn once in the bottom left panel
    add_reduction_arrows(
        p, reduction = reduction, labels = dims$labels, arrow_args = arrow_args,
        panel = if (rd$split_plot) levels(df$panel)[layout$bottom_left]
    )
}


#' @title Save a fixed-panel plot sized to fit its contents
#'
#' @description Saves a plot whose panels have a fixed absolute size (e.g. from
#'              \code{\link{DimPlot_bioToolkit}} or
#'              \code{\link{FeaturePlot_bioToolkit}}), including patchwork
#'              combinations of such plots. The width and height of the file are
#'              computed from the panels, strips, legends and margins, measured
#'              with the same kind of device the file is written with so that
#'              font metrics match.
#'
#' @param plot A ggplot or patchwork object with fixed panel sizes
#'
#' @param filename File to save to; the format is inferred from the extension
#'
#' @param dpi Resolution for raster formats, default is 1000
#'
#' @param bg Background colour. If NULL (default), taken from the plot
#'
#' @param ... Additional arguments passed to \code{\link[ggplot2]{ggsave}}
#'
#' @return The width and height of the saved file in inches, invisibly
#'
#' @export
#'
save_dimplot <- function(plot, filename, dpi = 1000, bg = NULL, ...) {
    # Build the gtable once, on the same kind of device it will be saved to, and
    # draw that exact gtable so measured and saved layouts cannot differ
    raster_ext <- grepl("\\.(png|jpe?g|tiff?|bmp)$", filename, ignore.case = TRUE)
    old_dev <- grDevices::dev.cur()
    if (raster_ext && requireNamespace("ragg", quietly = TRUE)) {
        tmp <- tempfile(fileext = ".png")
        ragg::agg_png(tmp, res = dpi)
        on.exit(unlink(tmp), add = TRUE)
    } else {
        grDevices::pdf(NULL)
    }
    measure_dev <- grDevices::dev.cur()
    is_patchwork <- inherits(plot, "patchwork")
    gt <- if (is_patchwork) patchwork::patchworkGrob(fix_nested_patchworks(plot)) else ggplotGrob(plot)
    dims <- fixed_gtable_size(gt, patchwork = is_patchwork)
    grDevices::dev.off(measure_dev)
    if (old_dev > 1) grDevices::dev.set(old_dev)

    if (is.null(bg)) {
        background <- gt$grobs[[match("background", gt$layout$name)]]
        bg <- if (is.null(background$gp$fill)) "transparent" else background$gp$fill
    }

    ggsave(filename, gt, width = dims[["width"]], height = dims[["height"]], units = "in",
           dpi = dpi, bg = bg, limitsize = FALSE, ...)
    invisible(dims)
}


#' @rdname save_dimplot
#'
#' @description `save_umap()` is a deprecated alias of `save_dimplot()`, which
#'              is named for the fact that reductions other than UMAP are
#'              plotted just as well.
#'
#' @export
#'
save_umap <- function(plot, filename, dpi = 1000, bg = NULL, ...) {
    warning("`save_umap()` is deprecated, use `save_dimplot()` instead", call. = FALSE)
    save_dimplot(plot = plot, filename = filename, dpi = dpi, bg = bg, ...)
}



