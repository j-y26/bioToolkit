# Publication-quality feature plots on a reduction for Seurat objects

#' @title Publication-quality feature plot
#'
#' @description Colours cells on a 2D reduction (UMAP, tSNE, PCA, ...) by the
#'              values of one or more features (genes, metadata columns, module
#'              scores, ...). Like \code{\link{DimPlot_bioToolkit}}, every panel
#'              has a fixed absolute size set by `panel.width` and
#'              `panel.height`, so that panels are identical across figures; use
#'              \code{\link{save_plot}} to save the plot, or a patchwork
#'              combination of plots, sized to fit. Text sizes follow
#'              \code{\link{theme_bioToolkit}}, and the values
#'              behind each panel can be written out with
#'              \code{\link{export_source_data}}.
#'
#'              Cells go through three steps:
#'              \enumerate{
#'                \item `cells` selects the cells that are plotted at all; the
#'                      others are left out entirely.
#'                \item `split.by` divides the plotted cells into panels.
#'                \item Within each panel, cells that are in `cells.highlight`
#'                      and have a value above `na.cutoff` are coloured by
#'                      their value. All other plotted cells are drawn beneath
#'                      them in `na.value` grey, which is never shown in the
#'                      legend and does not affect the colour scale.
#'              }
#'
#'              Values can be shown per cell (`geom = "point"`), averaged into
#'              hexagonal tiles (`geom = "hex"`), or smoothed over neighbouring
#'              cells (`geom = "density"`), which recovers the signal of sparse
#'              genes.
#'
#' @param object A Seurat object
#'
#' @param features Features to plot: genes, metadata columns or anything else
#'                 \code{Seurat::FetchData()} accepts. Several features are
#'                 drawn as one panel each, and only one feature can be given
#'                 with `split.by`
#'
#' @param cells Cells to plot, as names or indices (e.g. from
#'              \code{Seurat::WhichCells()}). Other cells are left out of the
#'              plot entirely. If NULL (default), all cells are plotted
#'
#' @param cells.highlight Cells to colour, as names or indices. Other plotted
#'                        cells are drawn in `na.value`. If NULL (default), all
#'                        plotted cells can be coloured
#'
#' @param full.extent Whether the axes span all cells of the object rather than
#'                    only the plotted ones, so that a plot of a subset of cells
#'                    lines up with a plot of the whole object. Default is FALSE
#'
#' @param split.by Metadata column to split the plot by, one panel per level.
#'                 If NULL (default), a single panel per feature is drawn
#'
#' @param keep.scale How colour scales are shared between panels, as in
#'                   \code{Seurat::FeaturePlot()}: "feature" (default) gives all
#'                   panels of a feature one scale and one legend; "all" gives
#'                   every panel of every feature one scale; NULL gives each
#'                   panel its own scale and legend
#'
#' @param reduction Name of the reduction to plot, default is "umap"
#'
#' @param dims The two dimensions of the reduction to plot, default is 1 and 2
#'
#' @param assay Assay that genes are read from. If NULL (default), the
#'              default assay of the object (\code{SeuratObject::DefaultAssay()}),
#'              as in Seurat. The assay and layer used are reported in a
#'              message. Genes are fetched from this assay only: a gene that
#'              is only in another assay is an error rather than being read
#'              from it silently. Objects holding both "RNA" and "SCT" assays
#'              give different values for the same gene, so set it explicitly
#'              for reproducible figures
#'
#' @param slot Layer to fetch expression from, default is "data"
#'
#' @param geom How values are drawn: "point" (default) colours each cell; "hex"
#'             colours hexagonal tiles by a summary of the cells in them (see
#'             `hex.fun`); "density" colours each cell by values smoothed over
#'             neighbouring cells (see `density.method`)
#'
#' @param cols The colours of the scale. If NULL (default), "plasma_r" or,
#'             for a diverging scale, "blue-magenta". Otherwise two or more
#'             colours from low to high values, e.g. `c("grey90", "darkred")`
#'             (low, high) or, for a diverging scale, `c("#2166AC", "white",
#'             "#B2182B")` (low, centre, high); or the name of a sequential or
#'             diverging palette (see \code{\link{show_palettes}}), such as
#'             "viridis", "magma" or "tol-burd", or of an RColorBrewer
#'             palette
#'
#' @param na.value Colour of cells (or tiles) that are not coloured. If NULL
#'                 (default), "grey88", or "grey70" for a diverging scale so
#'                 that they stand out from its light centre
#'
#' @param na.cutoff Values at or below this are drawn in `na.value` rather than
#'                  on the colour scale. "auto" (default) uses 0 for features
#'                  with no negative values, so that cells with no expression
#'                  are told apart from cells with low expression, and no
#'                  cutoff otherwise, and for `geom = "density"`. NULL plots all
#'                  values on the colour scale
#'
#' @param min.cutoff,max.cutoff Lower and upper ends of the colour scale, one
#'                              per feature or a single value for all. Values
#'                              beyond them get the colour of the end. A number
#'                              sets the end exactly, so that separate plots
#'                              share colours; a quantile such as "q1" or "q99"
#'                              is computed from the coloured values; NA
#'                              (default) uses the smallest or largest value
#'
#' @param diverging Whether to use a diverging scale centred on `center`:
#'                  "auto" (default) does so when the coloured values span
#'                  `center`, or TRUE / FALSE
#'
#' @param center Centre of a diverging scale, default is 0. The two sides of
#'               the scale extend equally far in colour from the centre, so
#'               that a side with a smaller range does not reach its darkest
#'               colour
#'
#' @param order Whether to draw cells with higher values on top, default is
#'              TRUE; on a diverging scale, cells further from the centre are
#'              drawn on top. If FALSE, cells are drawn in random order
#'
#' @param pt.size Point size, default is 0.3
#'
#' @param alpha Point (or tile) transparency, default is 1
#'
#' @param outline Whether to draw an outline around the cell clouds, default is
#'                FALSE. Not used for `geom = "hex"`
#'
#' @param outline.color Colour of the outline, default is "black"
#'
#' @param outline.size Thickness of the outline (in point size units), default
#'                     is 0.6
#'
#' @param label Whether to label groups at their median position, default is
#'              FALSE. Groups are the active identities unless `label.by` is
#'              given
#'
#' @param label.by Metadata column whose groups are labelled when
#'                 `label = TRUE`. If NULL (default), the active identities
#'
#' @param label.size Size of the group labels in mm, default is 3.5
#'
#' @param label.color Colour of the group labels, default is "black"
#'
#' @param hex.bins Number of hexagons across the width of a panel, default is
#'                 60. Hexagons are regular for any panel shape
#'
#' @param hex.fun Summary of a tile: "mean" (default) is the mean value of its
#'                coloured-eligible cells, counting those at or below
#'                `na.cutoff`; "pct" is the percentage of them above
#'                `na.cutoff` (or 0 when there is none). Tiles with no such
#'                cell, or a summary at or below the cutoff, are drawn in
#'                `na.value`
#'
#' @param hex.min.cells Tiles with fewer plotted cells are not drawn, default
#'                      is 3
#'
#' @param density.method For `geom = "density"`: "smooth" (default) is a
#'                       kernel-weighted mean of the values of neighbouring
#'                       cells, which keeps the units of the feature and is
#'                       comparable between panels; "wkde" is a kernel density
#'                       of cells weighted by their values, as in the Nebulosa
#'                       package, scaled to a maximum of 1 per feature. It is
#'                       highest where many cells have high values, and needs
#'                       non-negative values
#'
#' @param density.adjust Multiplier of the smoothing bandwidth, default is 1.
#'                       Smaller values smooth less
#'
#' @param legend.position Position of the legend, one of "right" (default),
#'                        "left", "top", "bottom", "inside" or "none"
#'
#' @param legend.direction Direction of the colour bar, "vertical" or
#'                         "horizontal". If NULL (default), horizontal for top
#'                         and bottom legends and vertical otherwise
#'
#' @param legend.title Title of the colour bar: "auto" (default) for the
#'                     feature name (or "Expression" / "Value" for a scale
#'                     shared by several features); NULL for no title; or a
#'                     string, or one string per feature (optionally named by
#'                     feature)
#'
#' @param legend.title.face Font face of the legend title. "auto" (default)
#'                          sets it in bold, gene names in bold italics
#'
#' @param legend.length Length of the colour bar in inches, default is 1.5
#'
#' @param legend.width Thickness of the colour bar in inches, default is 0.15
#'
#' @param legend.text.size Font size of the legend text and title in points.
#'                         If NULL (default), the theme's base size (11)
#'
#' @param ncol Number of panel columns. If NULL (default), split panels are
#'             laid out in a single row, and several features in up to four
#'             columns as in \code{Seurat::FeaturePlot()}
#'
#' @param by.row Whether panels fill the grid row by row (default, TRUE) or
#'               column by column (FALSE)
#'
#' @param arrow.panels Panels that get the dimension arrows: "bottom.left"
#'                     (default; the first panel of the bottom row), "all",
#'                     "left" / "right" (the first / last panel of each row),
#'                     "top" / "bottom" (the top / lowest panel of each
#'                     column) or "none". With `combine = FALSE`, every plot
#'                     gets its own arrows unless "none"
#'
#' @param combine Whether to return a single plot (default, TRUE) or a list of
#'                one plot per panel, each with its own arrows, to be arranged
#'                by the caller
#'
#' @param panel.width Width of a single panel in inches, default is 4
#'
#' @param panel.height Height of a single panel in inches, default is
#'                     `panel.width`
#'
#' @param panel.spacing Gap between panels in inches, default is 0.1
#'
#' @param strip.text.size Font size of the panel labels in points. If NULL
#'                        (default), the theme's large size (12)
#'
#' @param strip.text.face Font face of the panel labels, default is "bold".
#'                        Gene names are additionally set in italics
#'
#' @param raster Whether to rasterize the points while keeping text, legend
#'               and arrows as vectors: NULL (default) rasterizes when there
#'               are more than `raster.threshold` cells, or TRUE / FALSE. Hex
#'               tiles are always drawn as vectors
#'
#' @param raster.threshold Number of cells above which `raster = NULL`
#'                         rasterizes, default is 10000
#'
#' @param raster.dpi Resolution of the rasterized points, default is 1000. Match
#'                   the `dpi` used when saving raster (e.g. png) output
#'
#' @param bg Background colour, default is "transparent", so that the figure
#'           can be placed on any background (gaps inside outlines and label
#'           halos are then white). Use e.g. "white" for an opaque background
#'
#' @param seed Random seed used for drawing order and label placement, default
#'             is 42. The global random number generator state is not changed
#'
#' @param ... Arguments of the dimension arrows, as in
#'            \code{\link{add_dimension_arrows}} (e.g. `arrow.length = 1.5`;
#'            the snake_case names used before 0.10.0 are deprecated)
#'
#' @return A ggplot object when there is a single colour scale over one
#'         feature; otherwise a patchwork object, whose panels keep their
#'         fixed size when combined further with patchwork. With
#'         `combine = FALSE`, a named list of ggplot objects
#'
#' @import ggplot2
#' @export
#'
FeaturePlot_bioToolkit <- function(
    object,
    features,
    cells = NULL,
    cells.highlight = NULL,
    full.extent = FALSE,
    split.by = NULL,
    keep.scale = "feature",
    reduction = "umap",
    dims = c(1, 2),
    assay = NULL,
    slot = "data",
    geom = c("point", "hex", "density"),
    cols = NULL,
    na.value = NULL,
    na.cutoff = "auto",
    min.cutoff = NA,
    max.cutoff = NA,
    diverging = "auto",
    center = 0,
    order = TRUE,
    pt.size = 0.3,
    alpha = 1,
    outline = FALSE,
    outline.color = "black",
    outline.size = 0.6,
    label = FALSE,
    label.by = NULL,
    label.size = 3.5,
    label.color = "black",
    hex.bins = 60,
    hex.fun = c("mean", "pct"),
    hex.min.cells = 3,
    density.method = c("smooth", "wkde"),
    density.adjust = 1,
    legend.position = c("right", "left", "top", "bottom", "inside", "none"),
    legend.direction = NULL,
    legend.title = "auto",
    legend.title.face = "auto",
    legend.length = 1.5,
    legend.width = 0.15,
    legend.text.size = NULL,
    ncol = NULL,
    by.row = TRUE,
    arrow.panels = c("bottom.left", "all", "left", "right", "top", "bottom", "none"),
    combine = TRUE,
    panel.width = 4,
    panel.height = panel.width,
    panel.spacing = 0.1,
    strip.text.size = NULL,
    strip.text.face = "bold",
    raster = NULL,
    raster.threshold = 1e4,
    raster.dpi = 1000,
    bg = "transparent",
    seed = 42,
    ...
) {

    # Preparation
    check_installed("Seurat", reason = "to fetch cell embeddings and features")
    check_fixed_panel_support()
    geom <- match.arg(geom)
    hex.fun <- match.arg(hex.fun)
    density.method <- match.arg(density.method)
    legend.position <- match.arg(legend.position)
    if (!is.null(keep.scale)) keep.scale <- match.arg(keep.scale, c("feature", "all"))
    if (!(identical(diverging, "auto") || isTRUE(diverging) || isFALSE(diverging))) {
        stop("`diverging` must be \"auto\", TRUE or FALSE", call. = FALSE)
    }
    if (is.null(legend.direction)) {
        legend.direction <- if (legend.position %in% c("top", "bottom")) "horizontal" else "vertical"
    }
    legend.direction <- match.arg(legend.direction, c("vertical", "horizontal"))
    if (geom == "hex") check_installed("hexbin", reason = "to bin cells into hexagons")
    arrow.panels <- match.arg(arrow.panels)
    arrow_args <- list(...)
    raster <- resolve_raster_arg(raster, arrow_args)
    arrow_args$rasterize <- NULL
    arrow_args <- resolve_arrow_args(arrow_args, "FeaturePlot_bioToolkit")
    dims <- resolve_dims(dims, reduction)
    # Text sizes not given are taken from theme_sizes(), as in theme_bioToolkit()
    sizes <- theme_sizes()
    if (is.null(legend.text.size)) legend.text.size <- sizes$base
    if (is.null(strip.text.size)) strip.text.size <- sizes$large

    features <- unique(features)
    n_features <- length(features)
    if (n_features == 0) stop("No `features` given", call. = FALSE)
    if (n_features > 1 && !is.null(split.by)) {
        stop("Only one feature can be plotted with `split.by`; plot one feature at a time, ",
             "e.g. with lapply(), and combine the plots with patchwork", call. = FALSE)
    }
    min.cutoff <- stats::setNames(rep_len(as.list(min.cutoff), n_features), features)
    max.cutoff <- stats::setNames(rep_len(as.list(max.cutoff), n_features), features)

    # Embedding of the plotted cells, with axis limits and split panels
    rd <- reduction_data(object, reduction, dims$dims, cells, split.by, full.extent)
    df <- rd$df
    df$highlight <- highlight_flags(object, rownames(df), cells.highlight)
    is_gene <- feature_is_gene(object, features)
    # The assay is resolved, and reported, only when genes are plotted:
    # metadata columns such as module scores do not come from an assay
    assay_info <- if (any(is_gene)) resolve_assay(object, assay, slot) else NULL
    values <- fetch_features(object, features, rownames(df), slot, assay_info)
    if (label) {
        group <- if (is.null(label.by)) {
            Seurat::Idents(object)[rownames(df)]
        } else {
            fetch_meta(object, label.by, rownames(df))
        }
        df$group <- if (is.factor(group)) group else factor(group)
    }

    # One panel per split level of the single feature, or per feature
    panels <- if (rd$split_plot) {
        lapply(seq_len(nlevels(df$panel)), function(i) list(
            feature = features, label = levels(df$panel)[i],
            rows = which(as.integer(df$panel) == i)
        ))
    } else {
        lapply(features, function(f) list(feature = f, label = f, rows = seq_len(nrow(df))))
    }
    for (i in seq_along(panels)) {
        panels[[i]]$scale <- if (is.null(keep.scale)) {
            paste0("panel", i)
        } else if (keep.scale == "all") {
            "all"
        } else {
            panels[[i]]$feature
        }
    }

    # Values per cell, or per tile, of every panel
    na_cutoffs <- lapply(features, function(f) resolve_na_cutoff(na.cutoff, values[[f]], geom))
    names(na_cutoffs) <- features
    bandwidth <- if (geom == "density") density_bandwidth(df$x, df$y, density.adjust)
    shape <- panel.height / panel.width
    for (i in seq_along(panels)) {
        pn <- panels[[i]]
        cut <- na_cutoffs[[pn$feature]]
        v <- values[[pn$feature]][pn$rows]
        cells_df <- df[pn$rows, , drop = FALSE]
        panels[[i]]$cells <- cells_df
        panels[[i]]$marks <- switch(
            geom,
            point = point_marks(cells_df, v, cut),
            density = density_marks(cells_df, v, cut, rd, bandwidth, density.method),
            hex = hex_marks(cells_df, v, cut, rd, hex.bins, shape, hex.fun, hex.min.cells)
        )
    }
    # The weighted density has no units of its own and is scaled to a maximum
    # of 1 per feature, over all of its panels so that they stay comparable
    if (geom == "density" && density.method == "wkde") {
        for (f in features) {
            idx <- which(vapply(panels, function(pn) pn$feature == f, logical(1)))
            top <- max(unlist(lapply(panels[idx], function(pn) pn$marks$value)), 0, na.rm = TRUE)
            if (top > 0) for (i in idx) panels[[i]]$marks$value <- panels[[i]]$marks$value / top
        }
    }

    # Colour scale of each group of panels sharing one
    suffix <- switch(
        geom,
        point = NULL,
        hex = if (hex.fun == "pct") "(% cells)",
        density = if (density.method == "smooth") "(smoothed)" else "(density)"
    )
    colour_scales <- list()
    for (key in unique(vapply(panels, `[[`, "", "scale"))) {
        members <- Filter(function(pn) pn$scale == key, panels)
        scale_features <- unique(vapply(members, `[[`, "", "feature"))
        colour_scales[[key]] <- resolve_colour_scale(
            members, scale_features, min.cutoff, max.cutoff, diverging, center,
            cols, na.value, key
        )
        colour_scales[[key]]$title <- legend_title_label(
            legend.title, legend.title.face, scale_features, features, is_gene, suffix
        )
    }

    # Draw. One ggplot when a single scale covers a single feature (facetted
    # when split); otherwise one ggplot per panel, combined with patchwork
    halo.color <- halo_colour(bg)
    build <- function(members, facet_ncol, facet_nrow, facet_dir, strip) {
        p <- feature_ggplot(
            members = members, colour_scale = colour_scales[[members[[1]]$scale]],
            strip = strip, strip.face = strip_face(strip.text.face, is_gene[[members[[1]]$feature]], rd$split_plot),
            facet_ncol = facet_ncol, facet_nrow = facet_nrow, facet_dir = facet_dir,
            geom = geom, order = order, center = center, pt.size = pt.size, alpha = alpha,
            outline = outline, outline.color = outline.color, outline.size = outline.size,
            label = label, label.size = label.size, label.color = label.color,
            legend.position = legend.position, legend.direction = legend.direction,
            legend.length = legend.length, legend.width = legend.width,
            legend.text.size = legend.text.size, legend.title.face = legend.title.face,
            panel.width = panel.width, panel.height = panel.height, panel.spacing = panel.spacing,
            strip.text.size = strip.text.size, rd = rd,
            raster = raster, raster.threshold = raster.threshold, raster.dpi = raster.dpi,
            bg = bg, halo.color = halo.color, seed = seed
        )
        # The source data of a gene names the assay and layer it was read from
        f <- members[[1]]$feature
        if (!is.null(assay_info) && is_gene[[f]]) {
            sd <- attr(p, "source_data")
            sd$assay <- assay_info$assay
            sd$layer <- assay_info$layer
            attr(p, "source_data") <- sd
        }
        p
    }
    arrows <- function(p, panel = NULL) {
        add_reduction_arrows(p, reduction = reduction, labels = dims$labels,
                             arrow_args = arrow_args, panel = panel)
    }

    default_ncol <- if (rd$split_plot) length(panels) else seurat_feature_ncol(n_features)
    layout <- panel_layout(length(panels), ncol, by.row, default_ncol)
    single_plot <- combine && n_features == 1 && length(colour_scales) == 1
    # Panels that get the arrows (`arrow.panels`), in the order they fill the
    # grid
    idx <- arrow_panel_index(length(panels), layout$ncol, layout$nrow, layout$by.row, arrow.panels)

    if (single_plot) {
        p <- build(panels, layout$ncol, layout$nrow, layout$dir, strip = rd$split_plot)
        labels <- vapply(panels, `[[`, "", "label")
        return(arrows(p, panel = if (arrow.panels == "none") character(0) else if (rd$split_plot) labels[idx]))
    }

    # A strip names the panel whenever there is more than one
    plots <- lapply(panels, function(pn) build(list(pn), 1L, 1L, "h", strip = length(panels) > 1))
    names(plots) <- vapply(panels, `[[`, "", "label")
    if (!combine) {
        if (arrow.panels == "none") return(plots)
        return(lapply(plots, arrows))
    }

    plots[idx] <- lapply(plots[idx], arrows)
    # Panels sharing a scale have identical legends, which patchwork collects
    # into one; panels with their own scale keep their own legend
    collect <- length(colour_scales) == 1
    # A single row (or column) is given only its number of rows (or columns),
    # so that `|` (or `/`) extends it with further plots rather than nesting it
    grid <- if (layout$nrow == 1) {
        list(nrow = 1)
    } else if (layout$ncol == 1) {
        list(ncol = 1)
    } else {
        list(ncol = layout$ncol, nrow = layout$nrow, byrow = by.row)
    }
    combined <- do.call(patchwork::wrap_plots, c(
        list(plots, guides = if (collect) "collect" else "keep"), grid
    )) & theme(legend.position = legend.position)
    # A patchwork draws its own background, from ggplot2's global theme unless
    # set, which would hide the panels' (transparent by default) background
    transparent_bg <- is.na(bg) || identical(bg, "transparent")
    combined + patchwork::plot_annotation(theme = theme(
        plot.background = element_rect(fill = if (transparent_bg) "transparent" else bg, colour = NA)
    ))
}


#' Draw the panels sharing one colour scale as a single ggplot
#'
#' @param members Panels to draw, each with its `cells`, `marks` and `label`
#'
#' @param colour_scale Resolved colour scale, from
#'                     \code{\link{resolve_colour_scale}}
#'
#' @param strip Whether to label the panels with a strip
#'
#' @param strip.face Font face of the strip text
#'
#' @param facet_ncol,facet_nrow,facet_dir Facet layout
#'
#' @param rd Output of \code{\link{reduction_data}}
#'
#' @param halo.color Colour of the gap inside outlines and of label halos
#'
#' @inheritParams FeaturePlot_bioToolkit
#'
#' @return a ggplot object
#'
#' @keywords internal
feature_ggplot <- function(members,
                           colour_scale,
                           strip,
                           strip.face,
                           facet_ncol,
                           facet_nrow,
                           facet_dir,
                           geom,
                           order,
                           center,
                           pt.size,
                           alpha,
                           outline,
                           outline.color,
                           outline.size,
                           label,
                           label.size,
                           label.color,
                           legend.position,
                           legend.direction,
                           legend.length,
                           legend.width,
                           legend.text.size,
                           legend.title.face,
                           panel.width,
                           panel.height,
                           panel.spacing,
                           strip.text.size,
                           rd,
                           raster,
                           raster.threshold,
                           raster.dpi,
                           bg,
                           halo.color,
                           seed) {
    lv <- vapply(members, `[[`, "", "label")
    with_panel <- function(d, i) {
        d$panel <- factor(rep(lv[i], nrow(d)), levels = lv)
        d
    }
    marks <- do.call(rbind, lapply(seq_along(members), function(i) with_panel(members[[i]]$marks, i)))
    cells <- do.call(rbind, lapply(seq_along(members), function(i) with_panel(members[[i]]$cells, i)))
    marks$coloured <- marks$coloured & !is.na(marks$value)
    na_value <- colour_scale$na.value

    p <- ggplot(cells, aes(x = x, y = y))
    hex <- geom == "hex"

    if (hex) {
        # Tiles are drawn as polygons with an outline of their own colour, which
        # closes the hairline seams anti-aliasing leaves between them
        grey <- marks[!marks$coloured, ]
        if (nrow(grey)) {
            p <- p + geom_polygon(
                aes(group = tile), data = grey, fill = na_value, colour = na_value,
                alpha = alpha, linewidth = 0.1
            )
        }
        p <- p + geom_polygon(
            aes(group = tile, fill = value, colour = value), data = marks[marks$coloured, ],
            alpha = alpha, linewidth = 0.1
        )
    } else {
        layer_fun <- points_layer_fun(raster, nrow(marks), raster.threshold, raster.dpi)
        if (outline) {
            p <- p + outline_layers(marks, pt.size, outline.size, outline.color, halo.color, layer_fun)
        }
        grey <- marks[!marks$coloured, ]
        coloured <- marks[marks$coloured, ]
        # Ties are broken at random, so that no panel or region is
        # systematically drawn on top
        coloured <- with_seed(seed, coloured[sample.int(nrow(coloured)), ])
        if (order) {
            key <- if (colour_scale$diverging) abs(coloured$value - center) else coloured$value
            coloured <- coloured[order(key), ]
        }
        if (nrow(grey)) {
            p <- p + layer_fun(geom_point(data = grey, colour = na_value, size = pt.size, alpha = alpha, stroke = 0))
        }
        p <- p + layer_fun(geom_point(
            aes(colour = value), data = coloured, size = pt.size, alpha = alpha, stroke = 0
        ))
    }

    vertical <- legend.direction == "vertical"
    p <- p +
        scale_colour_gradientn(
            colours = colour_scale$colours,
            limits = colour_scale$limits,
            oob = scales::oob_squish,
            rescaler = colour_scale$rescaler,
            # A horizontal bar has room for fewer labels side by side
            breaks = if (vertical) waiver() else fitted_breaks(legend.length, legend.text.size),
            aesthetics = if (hex) c("colour", "fill") else "colour",
            name = colour_scale$title,
            guide = guide_colourbar(
                direction = legend.direction,
                theme = theme(
                    legend.key.height = unit(if (vertical) legend.length else legend.width, "in"),
                    legend.key.width = unit(if (vertical) legend.width else legend.length, "in"),
                    legend.frame = element_rect(colour = "black", linewidth = 0.3),
                    legend.ticks = element_line(colour = "black", linewidth = 0.3),
                    legend.ticks.length = unit(0.03, "in"),
                    legend.title.position = "top"
                )
            )
        ) +
        extent_scales(rd) +
        embedding_theme(
            panel.width = panel.width, panel.height = panel.height,
            ncol = facet_ncol, nrow = facet_nrow, panel.spacing = panel.spacing,
            legend.position = legend.position, legend.text.size = legend.text.size,
            strip.text.size = strip.text.size, strip.text.face = strip.face, bg = bg
        ) +
        theme(
            legend.title = element_text(
                size = legend.text.size,
                face = if (identical(legend.title.face, "auto")) "bold" else legend.title.face,
                hjust = if (vertical) 0 else 0.5,
                margin = margin(b = 0.1, unit = "cm")
            )
        )

    if (strip) p <- p + facet_wrap(~ panel, ncol = facet_ncol, nrow = facet_nrow, dir = facet_dir)

    if (label) {
        p <- p + group_label_layer(
            cells, by_panel = strip, label.size = label.size,
            label.color = label.color, halo.color = halo.color, seed = seed
        )
    }
    # The values behind the panels, for export_source_data(): one row per cell
    # with the value it is coloured by, or one row per tile for hexagons
    attr(p, "source_data") <- feature_source_data(marks, hex)
    p
}


#' Values behind the panels of a feature plot
#'
#' @param marks Marks of all panels, with a panel column
#'
#' @param hex Whether the marks are hexagon vertices
#'
#' @return For cells, a data frame with cell, panel, x, y, value and coloured
#'         columns; for hexagons, one row per tile with its centre
#'
#' @keywords internal
feature_source_data <- function(marks,
                                hex) {
    if (!hex) return(marks[, c("cell", "panel", "x", "y", "value", "coloured")])
    key <- paste(marks$panel, marks$tile)
    first <- !duplicated(key)
    data.frame(
        panel = marks$panel[first], tile = marks$tile[first],
        x = stats::ave(marks$x, key)[first], y = stats::ave(marks$y, key)[first],
        value = marks$value[first], coloured = marks$coloured[first]
    )
}


#' Colour bar breaks whose labels fit side by side along a bar
#'
#' Breaks are chosen as usual, then every other one is dropped until their
#' labels, estimated from their number of characters, fit along the bar.
#'
#' @param length Length of the bar in inches
#'
#' @param text_size Font size of the labels
#'
#' @return A function of the limits returning the breaks
#'
#' @keywords internal
fitted_breaks <- function(length,
                          text_size) {
    function(limits) {
        breaks <- scales::breaks_extended(5)(limits)
        breaks <- breaks[breaks >= limits[1] - 1e-9 & breaks <= limits[2] + 1e-9]
        if (length(breaks) < 3) return(breaks)
        # Width of one label plus a gap, from an average character width of
        # about 0.55 em
        slot <- (max(nchar(format(breaks, trim = TRUE))) + 1.5) * text_size * 0.55 / 72
        max_n <- max(2, floor(length / slot))
        while (length(breaks) > max_n) breaks <- breaks[seq(1, length(breaks), by = 2)]
        breaks
    }
}


#' Number of columns Seurat uses for several features
#'
#' @param n Number of features
#'
#' @keywords internal
seurat_feature_ncol <- function(n) {
    if (n == 1) 1L else if (n > 9) 4L else if (n > 6) 3L else 2L
}


#' Font face of panel strips, in italics for genes
#'
#' @param face Font face of the strip text
#'
#' @param is_gene Whether the feature is a gene
#'
#' @param split_plot Whether strips name split levels rather than features
#'
#' @keywords internal
strip_face <- function(face,
                       is_gene,
                       split_plot) {
    if (split_plot || !is_gene) return(face)
    switch(face, plain = "italic", bold = "bold.italic", face)
}


#' Resolve `na.cutoff` for one feature
#'
#' @param na.cutoff "auto", NULL or a number
#'
#' @param v Values of the feature over all plotted cells
#'
#' @param geom The geom used
#'
#' @return NULL or a number
#'
#' @keywords internal
resolve_na_cutoff <- function(na.cutoff,
                              v,
                              geom) {
    if (is.null(na.cutoff)) return(NULL)
    if (identical(na.cutoff, "auto")) {
        if (geom == "density" || any(v < 0, na.rm = TRUE)) return(NULL)
        return(0)
    }
    if (!is.numeric(na.cutoff) || length(na.cutoff) != 1) {
        stop("`na.cutoff` must be \"auto\", NULL or a single number", call. = FALSE)
    }
    na.cutoff
}


#' Whether values are above a cutoff (all are when there is none)
#'
#' @keywords internal
above_cutoff <- function(v,
                         cut) {
    if (is.null(cut)) rep(TRUE, length(v)) else !is.na(v) & v > cut
}


#' One mark per cell for `geom = "point"`
#'
#' @param cells Cells of a panel, with x, y and highlight columns
#'
#' @param v Values of the feature for these cells
#'
#' @param cut The resolved `na.cutoff`
#'
#' @return A data frame with cell, x, y, value and coloured columns
#'
#' @keywords internal
point_marks <- function(cells,
                        v,
                        cut) {
    data.frame(
        cell = rownames(cells), x = cells$x, y = cells$y, value = v,
        coloured = cells$highlight & !is.na(v) & above_cutoff(v, cut)
    )
}


#' One mark per cell, coloured by smoothed values, for `geom = "density"`
#'
#' Only highlighted cells with a value take part in the smoothing, and only they
#' are coloured.
#'
#' @inheritParams point_marks
#'
#' @param rd Output of \code{\link{reduction_data}}
#'
#' @param bandwidth Kernel standard deviations along x and y
#'
#' @param method "smooth" or "wkde"
#'
#' @keywords internal
density_marks <- function(cells,
                          v,
                          cut,
                          rd,
                          bandwidth,
                          method) {
    use <- cells$highlight & !is.na(v)
    if (method == "wkde" && any(v[use] < 0)) {
        stop("`density.method = \"wkde\"` needs non-negative values; use \"smooth\"", call. = FALSE)
    }
    value <- rep(NA_real_, nrow(cells))
    if (any(use)) {
        value[use] <- smooth_values(cells$x[use], cells$y[use], v[use], rd$xlim, rd$ylim, bandwidth, method)
    }
    data.frame(
        cell = rownames(cells), x = cells$x, y = cells$y, value = value,
        coloured = use & above_cutoff(value, cut)
    )
}


#' Smoothing bandwidth for `geom = "density"`
#'
#' Scott's rule in two dimensions, with a robust estimate of the spread so that
#' the many separate clusters of an embedding do not inflate it.
#'
#' @param x,y Coordinates of all plotted cells
#'
#' @param adjust Multiplier of the bandwidth
#'
#' @return Kernel standard deviations along x and y
#'
#' @keywords internal
density_bandwidth <- function(x,
                              y,
                              adjust) {
    n <- length(x)
    bw <- function(z) {
        s <- min(stats::sd(z), stats::IQR(z) / 1.349)
        if (!is.finite(s) || s <= 0) s <- stats::sd(z)
        if (!is.finite(s) || s <= 0) s <- 1
        s * n^(-1 / 6)
    }
    adjust * c(bw(x), bw(y))
}


#' Kernel-smoothed values at the cells' own positions
#'
#' The cells' values (and, for "smooth", their counts) are smoothed on a grid
#' (\code{kde_grid()}) and read back at every cell by bilinear interpolation
#' (\code{grid_interpolate()}).
#'
#' @param x,y Coordinates of the cells
#'
#' @param w Values of the cells
#'
#' @param xlim,ylim Axis limits, which the grid spans
#'
#' @param bandwidth Kernel standard deviations along x and y
#'
#' @param method "smooth" for the kernel-weighted mean of `w`; "wkde" for the
#'               kernel density weighted by `w`, per cell
#'
#' @param n_grid Number of grid points along each axis, default is 256
#'
#' @return A numeric vector, one value per cell
#'
#' @keywords internal
smooth_values <- function(x,
                          y,
                          w,
                          xlim,
                          ylim,
                          bandwidth,
                          method,
                          n_grid = 256L) {
    weighted <- grid_interpolate(kde_grid(x, y, w, xlim, ylim, bandwidth, n_grid), x, y, xlim, ylim)
    if (method == "wkde") return(weighted / length(w))
    # Interpolating the two sums before dividing keeps the ratio stable where
    # few cells are near a grid point
    counts <- kde_grid(x, y, rep(1, length(w)), xlim, ylim, bandwidth, n_grid)
    weighted / grid_interpolate(counts, x, y, xlim, ylim)
}


#' Gaussian kernel sum of weighted cells on a regular grid
#'
#' Cells are binned onto a regular grid over the axis limits (nearest grid
#' point), and the binned weights are convolved with a Gaussian kernel,
#' applied along x and along y as two matrix products. With unit weights, the
#' result is a kernel density estimate up to a constant factor. Shared by
#' `geom = "density"` of \code{\link{FeaturePlot_bioToolkit}} and by
#' \code{\link{DensityDiffPlot_bioToolkit}}.
#'
#' @param x,y Coordinates of the cells
#'
#' @param w Weights of the cells
#'
#' @param xlim,ylim Axis limits, which the grid spans
#'
#' @param bandwidth Kernel standard deviations along x and y
#'
#' @param n_grid Number of grid points along each axis
#'
#' @return An `n_grid` x `n_grid` matrix: rows along x, columns along y
#'
#' @keywords internal
kde_grid <- function(x,
                     y,
                     w,
                     xlim,
                     ylim,
                     bandwidth,
                     n_grid = 256L) {
    gx <- seq(xlim[1], xlim[2], length.out = n_grid)
    gy <- seq(ylim[1], ylim[2], length.out = n_grid)
    to_bin <- function(z, lim) {
        pmin(pmax(round((z - lim[1]) / diff(lim) * (n_grid - 1)) + 1, 1), n_grid)
    }
    bin <- to_bin(x, xlim) + (to_bin(y, ylim) - 1) * n_grid
    kx <- stats::dnorm(outer(gx, gx, "-"), sd = bandwidth[1])
    ky <- stats::dnorm(outer(gy, gy, "-"), sd = bandwidth[2])
    grid <- matrix(0, n_grid, n_grid)
    if (length(w) == 0) return(grid)
    sums <- rowsum(w, bin)
    grid[as.integer(rownames(sums))] <- sums
    kx %*% grid %*% ky
}


#' Values of a grid at given positions, by bilinear interpolation
#'
#' @param grid A square matrix from \code{kde_grid()}
#'
#' @param x,y Positions
#'
#' @param xlim,ylim Axis limits the grid spans
#'
#' @return A numeric vector along `x`
#'
#' @keywords internal
grid_interpolate <- function(grid,
                             x,
                             y,
                             xlim,
                             ylim) {
    n_grid <- nrow(grid)
    fx <- (x - xlim[1]) / diff(xlim) * (n_grid - 1) + 1
    fy <- (y - ylim[1]) / diff(ylim) * (n_grid - 1) + 1
    i <- pmin(pmax(floor(fx), 1), n_grid - 1)
    j <- pmin(pmax(floor(fy), 1), n_grid - 1)
    tx <- fx - i
    ty <- fy - j
    grid[cbind(i, j)] * (1 - tx) * (1 - ty) + grid[cbind(i + 1, j)] * tx * (1 - ty) +
        grid[cbind(i, j + 1)] * (1 - tx) * ty + grid[cbind(i + 1, j + 1)] * tx * ty
}


#' Hexagonal tiles for `geom = "hex"`
#'
#' The hexagon grid spans the axis limits and is shaped to the panel, so that
#' hexagons are regular on the page and identical across panels.
#'
#' @inheritParams point_marks
#'
#' @param rd Output of \code{\link{reduction_data}}
#'
#' @param bins Number of hexagons across the width of a panel
#'
#' @param shape Height to width ratio of a panel
#'
#' @param fun "mean" or "pct"
#'
#' @param min_cells Tiles with fewer plotted cells are not drawn
#'
#' @return A data frame of hexagon vertices, six rows per tile, with tile,
#'         x, y, value and coloured columns
#'
#' @keywords internal
hex_marks <- function(cells,
                      v,
                      cut,
                      rd,
                      bins,
                      shape,
                      fun,
                      min_cells) {
    hb <- hexbin::hexbin(cells$x, cells$y, xbins = bins, shape = shape,
                         xbnds = rd$xlim, ybnds = rd$ylim, IDs = TRUE)
    tile <- hb@cID
    centres <- hexbin::hcell2xy(hb)
    use <- cells$highlight & !is.na(v)
    idx <- match(tile[use], hb@cell)
    n_tiles <- length(hb@cell)
    n_use <- tabulate(idx, n_tiles)
    tile_sum <- function(z) vapply(split(z, factor(idx, levels = seq_len(n_tiles))), sum, numeric(1))
    summary <- if (fun == "mean") {
        tile_sum(v[use]) / n_use
    } else {
        positive <- above_cutoff(v[use], if (is.null(cut)) 0 else cut)
        100 * tile_sum(as.numeric(positive)) / n_use
    }
    summary[n_use == 0] <- NA
    coloured <- n_use > 0 & (if (fun == "mean") above_cutoff(summary, cut) else summary > 0)
    keep <- hb@count >= min_cells

    # Hexagon vertices, as in hexbin's own plotting
    sx <- bins / diff(rd$xlim)
    sy <- bins * shape / diff(rd$ylim)
    dx <- 0.5 / sx
    dy <- (1 / sqrt(3)) / (2 * sy)
    off_x <- c(dx, 0, -dx, -dx, 0, dx)
    off_y <- c(dy, 2 * dy, dy, -dy, -2 * dy, -dy)
    n <- sum(keep)
    data.frame(
        tile = rep(hb@cell[keep], each = 6),
        x = rep(centres$x[keep], each = 6) + rep(off_x, n),
        y = rep(centres$y[keep], each = 6) + rep(off_y, n),
        value = rep(unname(summary[keep]), each = 6),
        coloured = rep(coloured[keep], each = 6) %in% TRUE
    )
}


#' Resolve the colour scale shared by a group of panels
#'
#' @param members Panels sharing the scale
#'
#' @param scale_features Features of these panels
#'
#' @param min.cutoff,max.cutoff Lists of cutoffs named by feature
#'
#' @param key Name of the scale, for messages
#'
#' @inheritParams FeaturePlot_bioToolkit
#'
#' @return A list with the colour `limits`, `colours`, `rescaler`,
#'         `na.value` and whether it is `diverging`
#'
#' @keywords internal
resolve_colour_scale <- function(members,
                                 scale_features,
                                 min.cutoff,
                                 max.cutoff,
                                 diverging,
                                 center,
                                 cols,
                                 na.value,
                                 key) {
    feature_values <- lapply(scale_features, function(f) {
        unlist(lapply(Filter(function(pn) pn$feature == f, members), function(pn) {
            m <- pn$marks
            m$value[m$coloured & !is.na(m$value)]
        }))
    })
    names(feature_values) <- scale_features
    all_values <- unlist(feature_values)

    if (length(all_values) == 0) {
        warning("No cells to colour for ", paste(scale_features, collapse = ", "),
                "; all are drawn in `na.value`", call. = FALSE)
        limits <- c(0, 1)
    } else {
        ends <- vapply(scale_features, function(f) c(
            cutoff_value(min.cutoff[[f]], feature_values[[f]], min),
            cutoff_value(max.cutoff[[f]], feature_values[[f]], max)
        ), numeric(2))
        limits <- c(min(ends[1, ], na.rm = TRUE), max(ends[2, ], na.rm = TRUE))
    }

    is_diverging <- if (identical(diverging, "auto")) {
        length(all_values) > 0 && min(all_values) < center && max(all_values) > center
    } else {
        diverging
    }
    # Diverging: the centre colour at `center`, each side stretched to its
    # own end, so that the extreme of the shorter side is not faded
    # (rescale_diverging(); user decision 2026-10-14)
    rescaler <- if (is_diverging) {
        rescale_diverging(center)
    } else {
        function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) scales::rescale(x, to = to, from = from)
    }

    list(
        limits = limits,
        colours = resolve_continuous_palette(cols, is_diverging),
        rescaler = rescaler,
        na.value = if (is.null(na.value)) (if (is_diverging) "grey70" else "grey88") else na.value,
        diverging = is_diverging
    )
}


#' Resolve a cutoff given as a number, a quantile ("q##") or NA
#'
#' @param cutoff The cutoff
#'
#' @param v Values the quantile or default is computed from
#'
#' @param default Function giving the end when `cutoff` is NA (min or max)
#'
#' @return A number, NA when `v` is empty
#'
#' @keywords internal
cutoff_value <- function(cutoff,
                         v,
                         default) {
    if (length(v) == 0) return(NA_real_)
    if (is.null(cutoff) || (length(cutoff) == 1 && is.na(cutoff))) return(default(v))
    if (is.numeric(cutoff)) return(cutoff)
    if (is.character(cutoff) && grepl("^q[0-9.]+$", cutoff)) {
        q <- as.numeric(substring(cutoff, 2)) / 100
        if (q < 0 || q > 1) stop("Quantile cutoff '", cutoff, "' must be between q0 and q100", call. = FALSE)
        return(unname(stats::quantile(v, q)))
    }
    stop("Cutoffs must be numbers, quantiles such as \"q99\", or NA", call. = FALSE)
}


#' Legend title for a colour scale
#'
#' @param legend.title,legend.title.face As in
#'                                       \code{\link{FeaturePlot_bioToolkit}}
#'
#' @param scale_features Features sharing the scale
#'
#' @param features All features plotted
#'
#' @param is_gene Logical vector named by feature
#'
#' @param suffix Text appended to automatic titles, or NULL
#'
#' @return NULL, a string, or a plotmath expression with gene names in italics
#'
#' @keywords internal
legend_title_label <- function(legend.title,
                               legend.title.face,
                               scale_features,
                               features,
                               is_gene,
                               suffix) {
    if (is.null(legend.title)) return(NULL)
    if (!identical(legend.title, "auto")) {
        if (length(legend.title) == 1) return(legend.title)
        f <- scale_features[1]
        if (!is.null(names(legend.title)) && f %in% names(legend.title)) return(legend.title[[f]])
        return(legend.title[[match(f, features)]])
    }
    if (length(scale_features) > 1) {
        base <- if (all(is_gene[scale_features])) "Expression" else "Value"
        return(paste(c(base, suffix), collapse = " "))
    }
    f <- scale_features
    if (identical(legend.title.face, "auto") && is_gene[[f]]) {
        if (is.null(suffix)) return(bquote(bolditalic(.(f))))
        return(bquote(bolditalic(.(f)) ~ bold(.(suffix))))
    }
    paste(c(f, suffix), collapse = " ")
}
