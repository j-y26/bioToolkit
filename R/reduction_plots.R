# Functions to plot better reduction plots for Seurat objects

#' @title Publication-quality dimensionality reduction plot
#'
#' @description Plots cells on a 2D reduction (UMAP, tSNE, PCA, ...) coloured by
#'              a grouping variable. The panel (the scatter itself) has a fixed
#'              absolute size set by `panel.width` and `panel.height`, so it is
#'              identical across figures regardless of legend size or position.
#'              Use \code{\link{save_umap}} to save it: the output file is sized
#'              to fit the panel, legend and margins exactly.
#'
#' @param object A Seurat object
#'
#' @param group.by Metadata column to colour cells by. If NULL, the active
#'                 identities are used
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
#'              FALSE
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
#' @param panel.width Width of the panel in inches, default is 4
#'
#' @param panel.height Height of the panel in inches, default is `panel.width`
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
    panel.width = 4,
    panel.height = panel.width,
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
        Seurat::FetchData(object, vars = group.by, cells = rownames(emb))[[1]]
    }

    # Construct a data frame with embedding coordinates and group information
    df <- data.frame(x = emb[, 1], y = emb[, 2], group = group)
    if (!is.factor(df$group)) df$group <- factor(df$group)

    # Resolve colours up front so outlines can reuse them
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
            # Fixed panel size: the scatter itself is identical across plots
            # regardless of legend size or position
            panel.widths = unit(panel.width, "in"),
            panel.heights = unit(panel.height, "in"),
            plot.background = element_rect(fill = if (transparent_bg) "transparent" else bg, colour = NA),
            # Extra bottom/left margin makes room for the arrow labels
            plot.margin = margin(0.2, 0.2, 0.6, 0.6, "cm")
        )

    if (label) {
        check_installed("ggrepel", reason = "to label groups")
        centers <- stats::aggregate(cbind(x, y) ~ group, data = df, FUN = stats::median)
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

    # Arrows in place of axes
    add_dimension_arrows(p, reduction = reduction, ...)
}


#' @title Save a fixed-panel plot sized to fit its contents
#'
#' @description Saves a plot whose panel has a fixed absolute size (e.g. from
#'              \code{\link{DimPlot_bioToolkit}}). The width and height of the
#'              file are computed from the panel, legend and margins, measured
#'              with the same kind of device the file is written with so that
#'              font metrics match.
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
save_umap <- function(plot, filename, dpi = 1000, bg = NULL, ...) {
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
