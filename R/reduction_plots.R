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
#' @param cols Colours for the groups, either named by group or in the order of
#'             the group levels. If NULL, the default ggplot2 hue palette is used
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
#' @param rasterize Whether to rasterize the points while keeping text, legend
#'                  and arrows as vectors: "auto" (default) rasterizes when there
#'                  are more than `raster.threshold` cells, or TRUE / FALSE
#'
#' @param raster.threshold Number of cells above which `rasterize = "auto"`
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
#' @param ... Additional arguments passed to \code{\link{add_dimension_arrows}}
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
    cols = NULL,
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
    panel.width = 4,
    panel.height = panel.width,
    panel.spacing = 0.1,
    strip.text.size = 12,
    strip.text.face = "bold",
    rasterize = "auto",
    raster.threshold = 1e4,
    raster.dpi = 1000,
    bg = "white",
    seed = 42,
    ...
) {

    # Preparation
    check_installed("Seurat", reason = "to fetch cell embeddings and metadata")
    if (utils::packageVersion("ggplot2") < "4.0.0") {
        stop("ggplot2 >= 4.0.0 is required for fixed panel sizes", call. = FALSE)
    }
    legend.position <- match.arg(legend.position)
    if (!(identical(rasterize, "auto") || isTRUE(rasterize) || isFALSE(rasterize))) {
        stop("`rasterize` must be \"auto\", TRUE or FALSE", call. = FALSE)
    }

    # Fetch the embedding for the specified reduction
    emb <- Seurat::Embeddings(object, reduction = reduction)
    if (ncol(emb) < 2) stop("Reduction '", reduction, "' has fewer than 2 dimensions", call. = FALSE)
    
    # Fetch groups for the embedding's own cells so rows are guaranteed to align
    group <- if (is.null(group.by)) {
        Seurat::Idents(object)[rownames(emb)]
    } else {
        fetch_meta(object, group.by, rownames(emb))
    }

    # Construct a data frame with embedding coordinates and group information
    df <- data.frame(x = emb[, 1], y = emb[, 2], group = group)
    if (!is.factor(df$group)) df$group <- factor(df$group)

    # One panel per level of the splitting variable. Levels are dropped to what
    # is actually present so that no empty panel is drawn and the layout below
    # matches the panels ggplot2 creates
    split_plot <- !is.null(split.by)
    if (split_plot) {
        if (length(split.by) != 1) stop("`split.by` must be a single metadata column", call. = FALSE)
        split <- fetch_meta(object, split.by, rownames(emb))
        if (!is.factor(split)) split <- factor(split)
        split <- droplevels(split)
        # Cells with no value get their own panel rather than being dropped
        if (anyNA(split)) split <- addNA(split)
        df$split <- split
        n_panels <- nlevels(df$split)
        panel_ncol <- if (is.null(ncol)) n_panels else min(max(as.integer(ncol), 1L), n_panels)
        panel_nrow <- ceiling(n_panels / panel_ncol)
    } else {
        panel_ncol <- 1L
        panel_nrow <- 1L
    }

    # Resolve colours up front so outlines can reuse them. Levels are kept for
    # the whole data set, not per panel, so that every panel of a split plot is
    # coloured the same way and shares one complete legend
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
    na.color <- "grey50"

    # Shuffle and repel with a local seed; the caller's RNG state is restored
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
        on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    } else {
        on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE)
    }
    set.seed(seed)
    if (shuffle) df <- df[sample.int(nrow(df)), ]

    # Rasterize only the point layers; text, legend and arrows stay vector
    do_raster <- if (identical(rasterize, "auto")) nrow(df) > raster.threshold else rasterize
    if (do_raster) check_installed("ggrastr", reason = "to rasterize points")
    points_layer <- function(layer) {
        if (do_raster) ggrastr::rasterise(layer, dpi = raster.dpi) else layer
    }

    # Background transparency
    transparent_bg <- is.na(bg) || identical(bg, "transparent")

    # Colour of the gap inside the outline and of the label halo
    halo.color <- if (transparent_bg) "white" else bg

    p <- ggplot(df, aes(x = x, y = y))

    # Outline: a coloured layer then a background-coloured layer slightly larger
    # than the points, drawn beneath all cells, gives the cell clouds a crisp border
    if (outline) {
        if (identical(outline.color, "darken")) {
            check_installed("colorspace", reason = "to darken outline colours")
            dark <- colorspace::darken(c(cols[levels(df$group)], na.color), amount = outline.darken)
            df$outline <- unname(dark[as.integer(df$group)])
            df$outline[is.na(df$group)] <- dark[length(dark)]
            outline_layer <- geom_point(
                aes(colour = I(outline)), data = df,
                size = pt.size + outline.size, stroke = 0
            )
        } else {
            outline_layer <- geom_point(size = pt.size + outline.size, colour = outline.color, stroke = 0)
        }
        p <- p +
            points_layer(outline_layer) +
            points_layer(geom_point(size = pt.size + outline.size / 2, colour = halo.color, stroke = 0))
    }

    p <- p +
        points_layer(geom_point(aes(colour = group), size = pt.size, alpha = alpha, stroke = 0)) +
        # drop = FALSE keeps every group in the legend even when a panel, or the
        # whole plot, does not contain all of them
        scale_colour_manual(values = cols, na.value = na.color, drop = FALSE) +
        guides(colour = guide_legend(
            ncol = legend.ncol,
            override.aes = list(size = legend.pt.size, alpha = 1)
        )) +
        theme_void(base_size = legend.text.size) +
        theme(
            legend.position = legend.position,
            legend.title = element_blank(),
            legend.text = element_text(size = legend.text.size),
            legend.key.size = unit(0.45, "cm"),
            legend.key.spacing.y = unit(0.05, "cm"),
            legend.margin = margin(0, 0, 0, 0),
            # Gap between panel and legend on whichever side the legend is
            legend.box.spacing = unit(0.3, "cm"),
            legend.background = element_blank(),
            # Fixed panel size: each scatter is identical across plots regardless
            # of the number of panels and of the legend's size or position. Given
            # one unit per panel column / row, ggplot2 sizes each panel, rather
            # than the panel area as a whole
            panel.widths = unit(rep(panel.width, panel_ncol), "in"),
            panel.heights = unit(rep(panel.height, panel_nrow), "in"),
            panel.spacing = unit(panel.spacing, "in"),
            strip.text = element_text(
                size = strip.text.size, face = strip.text.face,
                margin = margin(b = 0.15, unit = "cm")
            ),
            strip.background = element_blank(),
            plot.background = element_rect(fill = if (transparent_bg) "transparent" else bg, colour = NA),
            # Extra bottom/left margin makes room for the arrow labels
            plot.margin = margin(0.2, 0.2, 0.6, 0.6, "cm")
        )

    if (split_plot) p <- p + facet_wrap(~ split, ncol = panel_ncol)

    if (label) {
        check_installed("ggrepel", reason = "to label groups")
        by <- if (split_plot) cbind(x, y) ~ group + split else cbind(x, y) ~ group
        centers <- stats::aggregate(by, data = df, FUN = stats::median)
        p <- p + ggrepel::geom_text_repel(
            data = centers,
            aes(label = group),
            size = label.size,
            fontface = "bold",
            colour = "black",
            bg.colour = halo.color,
            bg.r = 0.15,
            seed = seed,
            min.segment.length = Inf
        )
    }

    # Arrows in place of axes, drawn once in the bottom left panel
    if (split_plot) {
        lv <- levels(df$split)
        bottom_left <- factor(lv[(panel_nrow - 1) * panel_ncol + 1], levels = lv, exclude = NULL)
        apply_arrow_grob(
            p,
            dimension_arrow_grob(reduction = reduction, ...),
            data = data.frame(split = bottom_left)
        )
    } else {
        add_dimension_arrows(p, reduction = reduction, ...)
    }
}


#' @title Save a fixed-panel plot sized to fit its contents
#'
#' @description Saves a plot whose panels have a fixed absolute size (e.g. from
#'              \code{\link{DimPlot_bioToolkit}}). The width and height of the
#'              file are computed from the panels, strips, legend and margins,
#'              measured with the same kind of device the file is written with
#'              so that font metrics match.
#'
#' @param plot A ggplot object with a fixed panel size
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
    gt <- ggplotGrob(plot)
    dims <- fixed_gtable_size(gt)
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


#' Fetch one variable for a fixed set of cells
#'
#' `Seurat::FetchData()` silently drops cells whose value is `NA`, which would
#' misalign the result with the embedding it is meant to annotate. The values
#' are therefore indexed back onto `cells`, restoring the dropped ones as `NA`
#' so that they can be shown as their own group or panel.
#'
#' @param object A Seurat object
#'
#' @param var Name of the variable to fetch
#'
#' @param cells Cell names to fetch, in the order they are wanted
#'
#' @return A vector with one value per cell in `cells`
#'
#' @keywords internal
fetch_meta <- function(object, var, cells) {
    data <- withCallingHandlers(
        Seurat::FetchData(object, vars = var, cells = cells),
        warning = function(w) {
            if (grepl("missing data", conditionMessage(w))) invokeRestart("muffleWarning")
        }
    )
    data[[1]][match(cells, rownames(data))]
}


#' Absolute size of a gtable, including legends that overflow their cell
#'
#' A legend taller than a fixed-height panel (right/left) or wider than a
#' fixed-width panel (top/bottom) spills outside its cell and is not reflected
#' in the gtable's widths/heights. The legend is centred on its cell and the
#' gtable is centred on the device, so the device is enlarged to the legend size
#' plus twice the larger of the spaces before and after its cell, which keeps at
#' least the plot margin on both sides of the legend.
#'
#' @param gt A gtable with absolute widths and heights
#'
#' @return Named numeric vector with width and height in inches
#'
#' @keywords internal
fixed_gtable_size <- function(gt) {
    if (any(grid::unitType(gt$widths) == "null") || any(grid::unitType(gt$heights) == "null")) {
        stop("The plot has no fixed size; set theme(panel.widths, panel.heights)", call. = FALSE)
    }
    widths <- grid::convertWidth(gt$widths, "in", valueOnly = TRUE)
    heights <- grid::convertHeight(gt$heights, "in", valueOnly = TRUE)
    width <- sum(widths)
    height <- sum(heights)

    for (i in grep("^guide-box", gt$layout$name)) {
        legend <- gt$grobs[[i]]
        if (inherits(legend, "zeroGrob")) next
        pos <- gt$layout[i, ]
        legend_w <- grid::convertWidth(grid::grobWidth(legend), "in", valueOnly = TRUE)
        legend_h <- grid::convertHeight(grid::grobHeight(legend), "in", valueOnly = TRUE)
        before_w <- sum(widths[seq_len(pos$l - 1)])
        after_w <- sum(widths[-seq_len(pos$r)])
        before_h <- sum(heights[seq_len(pos$t - 1)])
        after_h <- sum(heights[-seq_len(pos$b)])
        if (legend_w > sum(widths[pos$l:pos$r]) + 1e-3) {
            width <- max(width, legend_w + 2 * max(before_w, after_w))
        }
        if (legend_h > sum(heights[pos$t:pos$b]) + 1e-3) {
            height <- max(height, legend_h + 2 * max(before_h, after_h))
        }
    }
    c(width = width, height = height)
}
