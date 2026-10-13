# Where on an embedding cells of one condition are enriched relative to
# another (roadmap item S4)
#
# Design
# ------
# `DensityDiffPlot_bioToolkit()` compares the density of cells of each
# condition on a 2D embedding, independently of how many cells each
# condition has:
#   1. One Gaussian kernel density per condition on a regular grid over the
#      axis limits (`kde_grid()`, shared with FeaturePlot_bioToolkit()'s
#      `geom = "density"`), with one bandwidth for all conditions (Scott's
#      rule on all plotted cells, `density_bandwidth()`), so that differences
#      do not come from different smoothing.
#   2. Each density is normalized to sum to 1, which removes the effect of
#      cell numbers: a condition with twice the cells but the same
#      distribution gives the same map.
#   3. Replicates (`sample.by`, optional): one normalized density per sample,
#      averaged within each condition, so that one large sample cannot drive
#      the map (every sample weighs the same, as in plot_composition_scatter()).
#   4. Measure: "difference", (p_compare - p_reference) / max(p_reference), so
#      that +0.5 means extra density of half the reference's peak; or
#      "log2ratio", log2((p_compare + e) / (p_reference + e)) with
#      e = 1% of the reference's peak, which keeps sparse regions from
#      producing extreme ratios.
#   5. Regions where both conditions are sparse are masked: grid points where
#      both densities are below the `mask.quantile` quantile of the mean
#      density at the cells' positions (5% of cells sit in sparser places).
#   6. Drawing: by default each cell is coloured by the value at its position
#      (`draw = "cells"`, the look of the feature plots), masked cells in
#      grey; or the grid itself as a raster over grey cells
#      (`draw = "raster"`). Diverging palette (default "blue-magenta"),
#      symmetric around 0, one scale for all panels.
#
# `split.by` is the condition column: one panel per compared level against
# `reference`, as DimPlot_bioToolkit()'s panels per level.
#
# Embedding distances are distorted (UMAP / tSNE), so the map shows where
# composition shifts, not how large they are in expression space; it is not
# a test. Differential abundance testing on samples (e.g. Milo, DA-seq) is the
# statistical counterpart.


#' @title Difference in cell density between conditions on an embedding
#'
#' @description Shows where on a 2D embedding (UMAP, tSNE, ...) the cells of a
#'              condition are enriched or depleted relative to a reference
#'              condition: the density of each condition is estimated with a
#'              Gaussian kernel, normalized so that cell numbers do not matter,
#'              and their difference is drawn on a diverging scale. One panel
#'              per compared condition, with a fixed size, as
#'              \code{\link{DimPlot_bioToolkit}}; save it with
#'              \code{\link{save_plot}}; the values behind it are returned by
#'              \code{\link{export_source_data}}.
#'
#' @param object A Seurat object
#'
#' @param split.by Metadata column of the conditions (e.g. genotype)
#'
#' @param reference Level of `split.by` the others are compared with. If NULL
#'                  (default), the first level
#'
#' @param compare Levels of `split.by` compared with `reference`, one panel
#'                each. If NULL (default), all other levels
#'
#' @param sample.by Metadata column of the samples. If given, one density is
#'                  computed per sample and the densities of a condition's
#'                  samples are averaged, so that every sample weighs the same
#'                  and one large sample cannot drive the map. If NULL
#'                  (default), the cells of a condition are pooled
#'
#' @param reduction Name of the reduction, default is "umap"
#'
#' @param dims Two dimensions of the reduction to plot, default is c(1, 2)
#'
#' @param cells Cells to use, as names or indices. If NULL (default), all
#'
#' @param full.extent Whether the axes span all cells of the object rather
#'                    than only the plotted ones, default is FALSE
#'
#' @param measure "difference" (default): the compared density minus the
#'                reference density, divided by the reference's peak density;
#'                or "log2ratio": the log2 ratio of the two densities, with 1\%
#'                of the reference's peak added to both
#'
#' @param draw "cells" (default): every cell coloured by the value at its
#'             position, masked cells in `na.value`; or "raster": the values
#'             on a grid, over the cells in `na.value`
#'
#' @param density.adjust Multiplier of the kernel bandwidth (Scott's rule on
#'                       all plotted cells), default is 1. Larger values
#'                       smooth more
#'
#' @param n.grid Number of grid points along each axis, default is 256
#'
#' @param mask.quantile Regions where both conditions are sparser than this
#'                      quantile of the density at the cells' positions are
#'                      masked (not coloured), default is 0.05
#'
#' @param cols Colours of the diverging scale: NULL (default) for
#'             "blue-magenta", a palette name, or colours from low to high,
#'             e.g. `c(low, mid, high)`
#'
#' @param limits The colour scale runs from -`limits` to `limits`. If NULL
#'               (default), the largest absolute value drawn
#'
#' @param na.value Colour of masked cells, default is "grey88"
#'
#' @param contour Whether to draw contour lines of the density of all plotted
#'                cells, default is FALSE
#'
#' @param pt.size Size of the cells, default is 0.3
#'
#' @param order Whether cells with the largest absolute values are drawn on
#'              top, default is TRUE (`draw = "cells"`)
#'
#' @param label Whether to label groups of cells at their median position,
#'              default is FALSE
#'
#' @param label.by Metadata column of the labelled groups. If NULL (default),
#'                 the active identities
#'
#' @param label.size,label.color Size and colour of the labels, default 3.5
#'                               and "black"
#'
#' @param legend.position Position of the legend: "right" (default), "left",
#'                        "top", "bottom" or "none"
#'
#' @param legend.title Title of the legend. If NULL (default), the measure and
#'                     the reference
#'
#' @param ncol Number of columns of panels. If NULL (default), all in one row
#'
#' @param panel.width,panel.height Size of each panel in inches, default 4;
#'                                 `panel.height` defaults to `panel.width`
#'
#' @param panel.spacing Gap between panels in inches, default is 0.1
#'
#' @param raster Whether to rasterize the cells: NULL (default) does so above
#'               `raster.threshold` cells per panel; or TRUE / FALSE
#'
#' @param raster.threshold Number of cells above which `raster = NULL`
#'                         rasterizes, default is 1e4
#'
#' @param raster.dpi Resolution of rasterized cells, default is 1000
#'
#' @param bg Background colour, default is "transparent"
#'
#' @param seed Random seed of the drawing order of ties and of label
#'             placement, default is 42
#'
#' @param ... Arguments passed to the dimension arrows, as in
#'            \code{\link{DimPlot_bioToolkit}}
#'
#' @details
#' Cells without a value of `split.by` (or `sample.by`) are left out, with a
#' message. All plotted cells set the bandwidth and the axis limits; cells of
#' conditions that are not compared are drawn but do not enter the densities.
#'
#' The map shows where the composition of the embedding shifts between
#' conditions, not whether the shift is significant: embedding distances are
#' distorted and cells are not replicates. Use a differential abundance test
#' on samples (e.g. Milo) for that. With few samples per condition, consider
#' `sample.by` so that one sample cannot dominate.
#'
#' The source data (\code{\link{export_source_data}}) hold, for
#' `draw = "cells"`, one row per cell and panel: `cell`, `panel`, `x`, `y`,
#' `condition`, `value`, the normalized densities `density_reference` and
#' `density_compare` at the cell, and whether it is `masked`; for
#' `draw = "raster"`, the same per grid point.
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' # Where are Mut cells enriched relative to WT?
#' DensityDiffPlot_bioToolkit(obj, split.by = "genotype", reference = "WT")
#'
#' # Every sample weighs the same; log2 ratio; group labels
#' p <- DensityDiffPlot_bioToolkit(obj, split.by = "genotype", reference = "WT",
#'                                 sample.by = "sample", measure = "log2ratio",
#'                                 label = TRUE, label.by = "cell_type")
#' save_plot(p, "density_difference.pdf")
#' }
#'
#' @seealso \code{\link{DimPlot_bioToolkit}}, \code{\link{FeaturePlot_bioToolkit}},
#'          \code{\link{plot_composition_scatter}}
#'
#' @import ggplot2
#' @export
#'
DensityDiffPlot_bioToolkit <- function(
    object,
    split.by,
    reference = NULL,
    compare = NULL,
    sample.by = NULL,
    reduction = "umap",
    dims = c(1, 2),
    cells = NULL,
    full.extent = FALSE,
    measure = c("difference", "log2ratio"),
    draw = c("cells", "raster"),
    density.adjust = 1,
    n.grid = 256,
    mask.quantile = 0.05,
    cols = NULL,
    limits = NULL,
    na.value = "grey88",
    contour = FALSE,
    pt.size = 0.3,
    order = TRUE,
    label = FALSE,
    label.by = NULL,
    label.size = 3.5,
    label.color = "black",
    legend.position = c("right", "left", "top", "bottom", "none"),
    legend.title = NULL,
    ncol = NULL,
    panel.width = 4,
    panel.height = panel.width,
    panel.spacing = 0.1,
    raster = NULL,
    raster.threshold = 1e4,
    raster.dpi = 1000,
    bg = "transparent",
    seed = 42,
    ...
) {
    check_installed("Seurat", reason = "to fetch cell embeddings")
    check_fixed_panel_support()
    measure <- match.arg(measure)
    draw <- match.arg(draw)
    legend.position <- match.arg(legend.position)
    if (!(is.numeric(mask.quantile) && length(mask.quantile) == 1 && mask.quantile >= 0 && mask.quantile < 1)) {
        stop("`mask.quantile` must be a single number in [0, 1)", call. = FALSE)
    }
    if (!is.null(limits) && !(is.numeric(limits) && length(limits) == 1 && limits > 0)) {
        stop("`limits` must be NULL or a single positive number", call. = FALSE)
    }
    if (length(split.by) != 1) stop("`split.by` must be a single metadata column", call. = FALSE)
    arrow_args <- list(...)
    dims <- resolve_dims(dims, reduction)
    sizes <- theme_sizes()

    # Cells, conditions and samples -------------------------------------------------------
    rd <- reduction_data(object, reduction, dims$dims, cells, NULL, full.extent)
    df <- rd$df
    df$condition <- fetch_meta(object, split.by, rownames(df))
    if (!is.null(sample.by)) df$sample <- fetch_meta(object, sample.by, rownames(df))
    keep <- !is.na(df$condition) & (if (is.null(sample.by)) TRUE else !is.na(df$sample))
    if (!all(keep)) {
        message(sum(!keep), " cell(s) without a value of `", paste(c(split.by, sample.by), collapse = "` or `"),
                "` are left out")
        df <- df[keep, , drop = FALSE]
    }
    df$condition <- droplevels(as.factor(df$condition))
    lv <- levels(df$condition)
    if (length(lv) < 2) stop("`split.by` has fewer than two levels in the cells used", call. = FALSE)
    if (is.null(reference)) reference <- lv[1]
    if (is.null(compare)) compare <- setdiff(lv, reference)
    bad <- setdiff(c(reference, compare), lv)
    if (length(bad)) {
        stop("Not level(s) of `", split.by, "`: ", paste(bad, collapse = ", "), " (levels: ",
             paste(lv, collapse = ", "), ")", call. = FALSE)
    }
    if (length(reference) != 1 || reference %in% compare) {
        stop("`reference` must be one level, not among `compare`", call. = FALSE)
    }

    # Densities ---------------------------------------------------------------------------------
    bandwidth <- density_bandwidth(df$x, df$y, density.adjust)
    dens <- lapply(stats::setNames(c(reference, compare), c(reference, compare)), function(cond) {
        condition_density(df[df$condition == cond, , drop = FALSE], rd, bandwidth, n.grid, !is.null(sample.by))
    })
    report <- vapply(names(dens), function(cond) {
        paste0(cond, " (", attr(dens[[cond]], "n_cells"), " cells",
               if (!is.null(sample.by)) paste0(", ", attr(dens[[cond]], "n_samples"), " samples"), ")")
    }, character(1))
    message("Density of ", paste(report, collapse = ", "), "; reference: ", reference)
    if (!is.null(sample.by) && any(vapply(dens, attr, 0, "n_samples") < 2)) {
        warning("A condition has a single sample: its density is that sample's", call. = FALSE)
    }
    ref <- dens[[reference]]
    peak <- max(ref)
    gx <- seq(rd$xlim[1], rd$xlim[2], length.out = n.grid)
    gy <- seq(rd$ylim[1], rd$ylim[2], length.out = n.grid)

    panels <- lapply(compare, function(cond) {
        cmp <- dens[[cond]]
        value <- if (measure == "difference") {
            (cmp - ref) / peak
        } else {
            eps <- 0.01 * peak
            log2((cmp + eps) / (ref + eps))
        }
        # Mask: both conditions sparser than the mask quantile of their mean
        # density at the positions of their cells
        both <- df$condition %in% c(reference, cond)
        at_cells <- grid_interpolate((ref + cmp) / 2, df$x[both], df$y[both], rd$xlim, rd$ylim)
        threshold <- stats::quantile(at_cells, mask.quantile, names = FALSE)
        masked <- pmax(ref, cmp) < threshold
        list(condition = cond, label = paste(cond, "vs", reference), value = value, masked = masked,
             cmp = cmp)
    })
    panel_lv <- vapply(panels, `[[`, "", "label")

    # Marks: cells coloured by the value at their position, or the grid ----------------------
    cells_all <- do.call(rbind, lapply(panels, function(pn) {
        v <- grid_interpolate(pn$value, df$x, df$y, rd$xlim, rd$ylim)
        m <- grid_interpolate(pn$masked * 1, df$x, df$y, rd$xlim, rd$ylim) > 0.5
        data.frame(
            cell = rownames(df), panel = factor(pn$label, levels = panel_lv), x = df$x, y = df$y,
            condition = df$condition, value = ifelse(m, NA_real_, v),
            density_reference = grid_interpolate(ref, df$x, df$y, rd$xlim, rd$ylim),
            density_compare = grid_interpolate(pn$cmp, df$x, df$y, rd$xlim, rd$ylim),
            masked = m, stringsAsFactors = FALSE
        )
    }))
    grid_all <- if (draw == "raster") {
        do.call(rbind, lapply(panels, function(pn) {
            g <- expand.grid(i = seq_len(n.grid), j = seq_len(n.grid))
            data.frame(
                panel = factor(pn$label, levels = panel_lv), x = gx[g$i], y = gy[g$j],
                value = ifelse(as.vector(pn$masked), NA_real_, as.vector(pn$value)),
                density_reference = as.vector(ref), density_compare = as.vector(pn$cmp),
                masked = as.vector(pn$masked)
            )
        }))
    }
    drawn <- if (draw == "cells") cells_all$value else grid_all$value
    lim <- if (is.null(limits)) max(abs(drawn), na.rm = TRUE) else limits
    if (!is.finite(lim) || lim == 0) lim <- 1

    # Plot ----------------------------------------------------------------------------------------
    n_panels <- length(panels)
    layout <- panel_layout(n_panels, ncol, TRUE, n_panels)
    layer_fun <- points_layer_fun(raster, nrow(df), raster.threshold, raster.dpi)
    grey <- cells_all[is.na(cells_all$value), , drop = FALSE]
    coloured <- cells_all[!is.na(cells_all$value), , drop = FALSE]
    if (draw == "raster") grey <- cells_all
    p <- ggplot(cells_all, aes(x = .data$x, y = .data$y))
    if (nrow(grey)) {
        p <- p + layer_fun(geom_point(data = grey, colour = na.value, size = pt.size, stroke = 0))
    }
    if (draw == "cells") {
        coloured <- with_seed(seed, coloured[sample.int(nrow(coloured)), , drop = FALSE])
        if (order) coloured <- coloured[order(abs(coloured$value)), , drop = FALSE]
        p <- p + layer_fun(geom_point(aes(colour = .data$value), data = coloured, size = pt.size, stroke = 0))
    } else {
        p <- p + geom_raster(aes(fill = .data$value), data = grid_all, interpolate = TRUE)
    }
    if (contour) {
        # Contours of the density of all plotted cells, at the levels that
        # 10%, 30%, 50%, 70% and 90% of cells lie above
        pooled <- kde_grid(df$x, df$y, rep(1, nrow(df)), rd$xlim, rd$ylim, bandwidth, n.grid)
        breaks <- stats::quantile(grid_interpolate(pooled, df$x, df$y, rd$xlim, rd$ylim),
                                  c(0.1, 0.3, 0.5, 0.7, 0.9), names = FALSE)
        g <- expand.grid(i = seq_len(n.grid), j = seq_len(n.grid))
        cg <- data.frame(x = gx[g$i], y = gy[g$j], z = as.vector(pooled))
        p <- p + geom_contour(aes(x = .data$x, y = .data$y, z = .data$z), data = cg, breaks = breaks,
                              inherit.aes = FALSE, colour = "grey35", linewidth = pt_to_linewidth(0.4))
    }

    if (is.null(legend.title)) {
        legend.title <- if (measure == "difference") {
            paste0("Density difference\n(vs ", reference, ")")
        } else {
            paste0("log2 density ratio\n(vs ", reference, ")")
        }
    }
    vertical <- !(legend.position %in% c("top", "bottom"))
    p <- p +
        scale_colour_gradientn(
            colours = resolve_continuous_palette(cols, diverging = TRUE),
            limits = c(-lim, lim), oob = scales::oob_squish, na.value = "transparent",
            aesthetics = c("colour", "fill"), name = legend.title,
            guide = guide_colourbar(
                direction = if (vertical) "vertical" else "horizontal",
                theme = theme(
                    legend.key.height = unit(if (vertical) 1.5 else 0.15, "in"),
                    legend.key.width = unit(if (vertical) 0.15 else 1.5, "in"),
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
            ncol = layout$ncol, nrow = layout$nrow, panel.spacing = panel.spacing,
            legend.position = legend.position, legend.text.size = sizes$base,
            strip.text.size = sizes$large, strip.text.face = "bold", bg = bg
        ) +
        theme(legend.title = element_text(size = sizes$base, hjust = if (vertical) 0 else 0.5,
                                          margin = margin(b = 0.1, unit = "cm"))) +
        facet_wrap(~ panel, ncol = layout$ncol, nrow = layout$nrow)

    if (label) {
        group <- if (is.null(label.by)) Seurat::Idents(object)[cells_all$cell] else
            fetch_meta(object, label.by, cells_all$cell)
        lab_df <- cells_all
        lab_df$group <- if (is.factor(group)) group else factor(group)
        p <- p + group_label_layer(lab_df, by_panel = TRUE, label.size = label.size,
                                   label.color = label.color, halo.color = halo_colour(bg), seed = seed)
    }
    p <- add_reduction_arrows(p, reduction = reduction, labels = dims$labels, arrow_args = arrow_args,
                              panel = panel_lv[layout$bottom_left])
    attr(p, "source_data") <- if (draw == "cells") cells_all else grid_all
    p
}


#' Normalized kernel density of one condition's cells on a grid
#'
#' @param d Cells of the condition (`x`, `y` and, per sample, `sample`)
#'
#' @param rd Result of \code{reduction_data()} (axis limits)
#'
#' @param bandwidth Kernel standard deviations
#'
#' @param n.grid Grid size
#'
#' @param per_sample Whether to average one density per sample
#'
#' @return A grid summing to 1, with attributes `n_cells` and `n_samples`
#'
#' @keywords internal
condition_density <- function(d,
                              rd,
                              bandwidth,
                              n.grid,
                              per_sample) {
    one <- function(x, y) {
        g <- kde_grid(x, y, rep(1, length(x)), rd$xlim, rd$ylim, bandwidth, n.grid)
        g / sum(g)
    }
    if (!per_sample) {
        out <- one(d$x, d$y)
        n_samples <- NA_integer_
    } else {
        by_sample <- split(seq_len(nrow(d)), droplevels(as.factor(d$sample)))
        out <- Reduce(`+`, lapply(by_sample, function(i) one(d$x[i], d$y[i]))) / length(by_sample)
        n_samples <- length(by_sample)
    }
    attr(out, "n_cells") <- nrow(d)
    attr(out, "n_samples") <- n_samples
    out
}
