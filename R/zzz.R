# Package infrastructure helpers
#
# `@import rlang` below keeps rlang's helpers (e.g. `%||%`, which base R only
# has from 4.4) available package-wide; it was on de_dot_plot(), removed in
# 0.7.0. The package's own check_installed() takes precedence over rlang's.

#' Base R functions used without a namespace prefix
#'
#' @importFrom grDevices hcl
#' @importFrom methods is
#' @importFrom stats coef cor.test lm median p.adjust
#' @importFrom utils write.csv
#' @keywords internal
#' @name bioToolkit-imports
NULL

# Column names used in dplyr / ggplot2 code (non-standard evaluation), declared
# so that R CMD check does not report them as undefined globals
utils::globalVariables(c(
    "Freq", "GeneName", "coef", "colour", "corr_FDR", "corr_p_value", "group",
    "p_value", "percent", "piValue", "ring", "tile", "value", "x", "y"
))


#' Check if a package is installed, with informative error message
#'
#' @param pkg Character string, the package name.
#' @param reason Optional character string explaining why the package is needed.
#'
#' @return Invisible TRUE if package is available; otherwise stops with an error.
#'
#' @import rlang
#'
#' @keywords internal
check_installed <- function(pkg,
                            reason = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- paste0("Package '", pkg, "' is required")
    if (!is.null(reason)) {
      msg <- paste0(msg, " ", reason)
    }
    # Determine install command
    bioc_pkgs <- c("edgeR", "GSEABase", "DOSE")
    if (pkg %in% bioc_pkgs) {
      msg <- paste0(
        msg,
        ".\nInstall it with: BiocManager::install('", pkg, "')"
      )
    } else {
      msg <- paste0(msg, ".\nInstall it with: install.packages('", pkg, "')")
    }
    stop(msg, call. = FALSE)
  }
  invisible(TRUE)
}


#' Accept deprecated argument names
#'
#' Older versions of a few functions (\code{format_sc_plots()},
#' \code{manipulate_sc_plot()}, \code{add_dimension_arrows()} and the DE table
#' helpers) used snake_case argument names. Since 0.10.0 they use the
#' package's dot.case names (see CONTRIBUTING.md). The old names are still
#' accepted: the caller collects them in `...` and passes them here, which
#' warns once (class "deprecatedWarning") and returns the values under their
#' new names, ready for \code{list2env(..., environment())} in the caller.
#'
#' @param dots `list(...)` of the caller
#'
#' @param renames Character vector of new argument names, named by the old
#'                names
#'
#' @param fn Name of the calling function, for messages
#'
#' @param given Names of the arguments the caller received (from
#'              \code{names(match.call())}), to refuse an old and a new name
#'              given together
#'
#' @return A named list of the values of the old arguments, named by the new
#'         names (empty when no old name was used)
#'
#' @keywords internal
renamed_args <- function(dots,
                         renames,
                         fn,
                         given = character()) {
    if (length(dots) == 0) return(list())
    old <- names(dots)
    # `...` only exists to catch old names, so anything else is a user error,
    # reported as R reports an unused argument
    if (is.null(old) || any(!nzchar(old))) {
        stop("unused unnamed argument in `", fn, "()`", call. = FALSE)
    }
    unknown <- setdiff(old, names(renames))
    if (length(unknown)) {
        stop("unused argument", if (length(unknown) > 1) "s", " in `", fn, "()`: ",
             paste0("`", unknown, "`", collapse = ", "), call. = FALSE)
    }
    new <- unname(renames[old])
    clash <- new %in% given
    if (any(clash)) {
        stop("In `", fn, "()`, give either `", old[clash][1], "` (deprecated) or `",
             new[clash][1], "`, not both", call. = FALSE)
    }
    msg <- paste0("In `", fn, "()`, ", paste0("`", old, "`", collapse = ", "),
                  if (length(old) > 1) " are" else " is",
                  " deprecated; use ", paste0("`", new, "`", collapse = ", "), " instead.")
    warning(warningCondition(msg, class = "deprecatedWarning"))
    stats::setNames(dots, new)
}


#' Resolve the arrow arguments passed through `...` of an embedding plot
#'
#' \code{DimPlot_bioToolkit()}, \code{FeaturePlot_bioToolkit()} and
#' \code{DensityDiffPlot_bioToolkit()} pass `...` to the dimension arrows. The
#' arguments take the dot.case names of \code{\link{add_dimension_arrows}}
#' (since 0.10.0); the old snake_case names are accepted with a deprecation
#' warning. `arrow.scale` (new in 1.0.0) has no old name. They are returned
#' under the snake_case names of the internal
#' \code{\link{dimension_arrow_grob}}.
#'
#' @param dots `list(...)` of the plot, without `rasterize`
#'
#' @param fn Name of the plot function, for messages
#'
#' @return A named list of arguments for \code{dimension_arrow_grob()}
#'
#' @keywords internal
resolve_arrow_args <- function(dots,
                               fn) {
    if (length(dots) == 0) return(list())
    old <- names(dots) %in% names(ARROW_ARG_RENAMES)
    new_dots <- c(dots[!old], renamed_args(dots[old], ARROW_ARG_RENAMES, fn, names(dots)[!old]))
    # Internal (snake_case) name of each accepted argument
    known <- c(ARROW_ARG_RENAMES, arrow_scale = "arrow.scale")
    unknown <- setdiff(names(new_dots), known)
    if (length(unknown) || is.null(names(new_dots))) {
        stop("unused argument in `", fn, "()`: ", paste0("`", unknown, "`", collapse = ", "),
             call. = FALSE)
    }
    if (!is.null(new_dots$arrow.scale)) check_arrow_scale(new_dots$arrow.scale)
    stats::setNames(new_dots, names(known)[match(names(new_dots), known)])
}


#' Check if an object is a ggplot
#'
#' @param x An object to check.
#'
#' @return Logical, TRUE if the object inherits from "gg".
#'
#' @keywords internal
is_ggplot <- function(x) {
  inherits(x, "gg")
}


#' Stop unless ggplot2 supports fixed panel sizes
#'
#' @keywords internal
check_fixed_panel_support <- function() {
    if (utils::packageVersion("ggplot2") < "4.0.0") {
        stop("ggplot2 >= 4.0.0 is required for fixed panel sizes", call. = FALSE)
    }
}


# Helpers shared by the reduction plots ---------------------------------------

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
fetch_meta <- function(object,
                       var,
                       cells) {
    data <- withCallingHandlers(
        Seurat::FetchData(object, vars = var, cells = cells),
        warning = function(w) {
            if (grepl("missing data", conditionMessage(w))) invokeRestart("muffleWarning")
        }
    )
    data[[1]][match(cells, rownames(data))]
}


#' A metadata column for a set of cells, as a factor with all of its levels
#'
#' Levels come from all cells of the object, not only from `cells`, so that
#' groups get the same colours in a plot of a subset of cells as in a plot of
#' all of them. Levels of a factor are kept; other values are sorted.
#'
#' @param object A Seurat object
#'
#' @param group.by Metadata column, or NULL for the active identities
#'
#' @param cells Cell names
#'
#' @return A factor along `cells`
#'
#' @keywords internal
fetch_group <- function(object,
                        group.by,
                        cells) {
    if (is.null(group.by)) return(Seurat::Idents(object)[cells])
    values <- fetch_meta(object, group.by, colnames(object))
    group <- if (is.factor(values)) values else factor(values)
    group[match(cells, colnames(object))]
}


#' Resolve the `raster` argument, accepting the deprecated `rasterize`
#'
#' @param raster NULL, TRUE or FALSE
#'
#' @param dots The list of `...` arguments, which may hold `rasterize`
#'
#' @return NULL, TRUE or FALSE
#'
#' @keywords internal
resolve_raster_arg <- function(raster,
                               dots) {
    if (!is.null(dots$rasterize)) {
        warning("`rasterize` is deprecated, use `raster` instead", call. = FALSE)
        raster <- if (identical(dots$rasterize, "auto")) NULL else dots$rasterize
    }
    if (!(is.null(raster) || isTRUE(raster) || isFALSE(raster))) {
        stop("`raster` must be NULL, TRUE or FALSE", call. = FALSE)
    }
    raster
}



#' Resolve the dimensions to plot and their arrow labels
#'
#' Before `dims` selected the dimensions, it was passed through `...` to the
#' arrows as their labels. A character `dims` is therefore still taken as the
#' arrow labels of the first two dimensions.
#'
#' @param dims Two dimension indices, or two arrow labels
#'
#' @param reduction Name of the reduction
#'
#' @return A list with the integer `dims` and their arrow `labels`
#'
#' @keywords internal
resolve_dims <- function(dims,
                         reduction) {
    if (is.character(dims)) return(list(dims = c(1L, 2L), labels = dims))
    if (length(dims) != 2 || anyNA(dims)) stop("`dims` must be two dimension indices", call. = FALSE)
    dims <- as.integer(dims)
    list(dims = dims, labels = reduction_dim_labels(reduction, dims))
}



#' Resolve cells given as names or indices
#'
#' @param cells NULL, a character vector of cell names, or integer indices
#'
#' @param all_cells All available cell names
#'
#' @param arg Name of the argument, for messages
#'
#' @return A character vector of cell names, a subset of `all_cells` in the
#'         order of `all_cells`
#'
#' @keywords internal
resolve_cells <- function(cells,
                          all_cells,
                          arg = "cells") {
    if (is.null(cells)) return(all_cells)
    if (is.list(cells)) cells <- unlist(cells, use.names = FALSE)
    if (is.numeric(cells)) cells <- all_cells[cells[cells >= 1 & cells <= length(all_cells)]]
    cells <- as.character(cells)
    found <- all_cells[all_cells %in% cells]
    n_missing <- length(setdiff(cells, all_cells))
    if (n_missing > 0) warning(n_missing, " cell(s) in `", arg, "` were not found", call. = FALSE)
    found
}



#' Flag the plotted cells that are highlighted
#'
#' @param object A Seurat object
#'
#' @param cells The plotted cells
#'
#' @param cells.highlight Cells to highlight as names or indices of the
#'                        object's cells, or NULL to highlight all
#'
#' @return A logical vector, one value per cell in `cells`
#'
#' @keywords internal
highlight_flags <- function(object,
                            cells,
                            cells.highlight) {
    if (is.null(cells.highlight)) return(rep(TRUE, length(cells)))
    highlight <- resolve_cells(cells.highlight, colnames(object), arg = "cells.highlight")
    flags <- cells %in% highlight
    if (!any(flags)) warning("None of the plotted cells are in `cells.highlight`", call. = FALSE)
    flags
}



#' Embedding, axis limits and split panels for the plotted cells
#'
#' @param object A Seurat object
#'
#' @param reduction Name of the reduction
#'
#' @param dims Two dimension indices
#'
#' @param cells Cells to plot, or NULL for all
#'
#' @param split.by Metadata column to split by, or NULL
#'
#' @param full.extent Whether the axis limits span all cells of the reduction
#'                    rather than only the plotted ones
#'
#' @return A list with `df` (x, y and, when split, panel per plotted cell,
#'         with cell names as row names), the axis limits `xlim` and `ylim`,
#'         `split_plot` and the number of panels `n_panels`
#'
#' @keywords internal
reduction_data <- function(object,
                           reduction,
                           dims,
                           cells,
                           split.by,
                           full.extent) {
    emb <- Seurat::Embeddings(object, reduction = reduction)
    if (ncol(emb) < max(dims)) {
        stop("Reduction '", reduction, "' has fewer than ", max(dims), " dimensions", call. = FALSE)
    }
    emb <- emb[, dims, drop = FALSE]
    keep <- resolve_cells(cells, rownames(emb))
    if (length(keep) == 0) stop("None of `cells` are in reduction '", reduction, "'", call. = FALSE)
    xy <- emb[keep, , drop = FALSE]
    df <- data.frame(x = xy[, 1], y = xy[, 2], row.names = keep)

    # Limits are fixed explicitly, with ggplot2's default 5% expansion, so that
    # they can span all cells and so that hexagons can be sized to them
    extent <- if (full.extent) emb else xy
    xlim <- scales::expand_range(range(extent[, 1]), mul = 0.05)
    ylim <- scales::expand_range(range(extent[, 2]), mul = 0.05)

    # One panel per level of the splitting variable. Levels are dropped to what
    # is actually present so that no empty panel is drawn and the layout
    # matches the panels ggplot2 creates
    split_plot <- !is.null(split.by)
    if (split_plot) {
        if (length(split.by) != 1) stop("`split.by` must be a single metadata column", call. = FALSE)
        split <- fetch_meta(object, split.by, keep)
        if (!is.factor(split)) split <- factor(split)
        split <- droplevels(split)
        # Cells with no value get their own panel rather than being dropped
        if (anyNA(split)) split <- addNA(split)
        df$panel <- split
    }

    list(
        df = df, xlim = xlim, ylim = ylim, split_plot = split_plot,
        n_panels = if (split_plot) nlevels(df$panel) else 1L
    )
}



#' Grid layout of a set of panels
#'
#' @param n Number of panels
#'
#' @param ncol Requested number of columns, or NULL for `default_ncol`
#'
#' @param by.row Whether panels fill the grid row by row (TRUE) or column by
#'               column (FALSE)
#'
#' @param default_ncol Number of columns used when `ncol` is NULL
#'
#' @return A list with `ncol`, `nrow`, the facet direction `dir` and the
#'         index of the bottom left panel, `bottom_left`
#'
#' @keywords internal
panel_layout <- function(n,
                         ncol,
                         by.row,
                         default_ncol) {
    ncol <- if (is.null(ncol)) default_ncol else as.integer(ncol)
    ncol <- min(max(ncol, 1L), n)
    nrow <- as.integer(ceiling(n / ncol))
    # Filling by column, the last columns may be left empty when n is not a
    # multiple of nrow, so the number of columns is that actually filled
    if (!by.row) ncol <- as.integer(ceiling(n / nrow))
    list(
        ncol = ncol, nrow = nrow, dir = if (by.row) "h" else "v",
        bottom_left = if (by.row) (nrow - 1L) * ncol + 1L else nrow,
        by.row = by.row
    )
}


# Values of `arrow.panels`: which panels of an embedding grid get the
# dimension arrows (revision of 2026-10-14)
ARROW_PANELS <- c("bottom.left", "all", "left", "right", "top", "bottom", "none")


#' Panels of a grid that get the dimension arrows
#'
#' Panels are numbered in the order they fill the grid (by row or by column).
#' The edges are those of the panels actually drawn: in a grid whose last row
#' is incomplete, a panel with no panel below it is on the "bottom" edge, and
#' the last panel of a row is on the "right" edge.
#'
#' @param n Number of panels
#'
#' @param ncol,nrow Size of the grid
#'
#' @param by.row Whether panels fill the grid row by row
#'
#' @param which One of `ARROW_PANELS`: "bottom.left" (the first panel of the
#'              bottom row), "all", "left" / "right" (first / last panel of
#'              each row), "top" / "bottom" (top / lowest panel of each
#'              column) or "none"
#'
#' @return Indices of the panels, in 1..n
#'
#' @keywords internal
arrow_panel_index <- function(n,
                              ncol,
                              nrow,
                              by.row,
                              which) {
    which <- match.arg(which, ARROW_PANELS)
    i <- seq_len(n)
    row <- if (by.row) (i - 1L) %/% ncol + 1L else (i - 1L) %% nrow + 1L
    col <- if (by.row) (i - 1L) %% ncol + 1L else (i - 1L) %/% nrow + 1L
    keep <- switch(
        which,
        all = rep(TRUE, n),
        none = rep(FALSE, n),
        bottom.left = row == max(row) & col == min(col[row == max(row)]),
        left = col == stats::ave(col, row, FUN = min),
        right = col == stats::ave(col, row, FUN = max),
        top = row == stats::ave(row, col, FUN = min),
        bottom = row == stats::ave(row, col, FUN = max)
    )
    i[keep]
}



#' Run code with a fixed random seed, restoring the caller's RNG state after
#'
#' @param seed Random seed
#'
#' @param code Code to run
#'
#' @return The value of `code`
#'
#' @keywords internal
with_seed <- function(seed,
                      code) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
        on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    } else {
        on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE)
    }
    set.seed(seed)
    code
}


#' A function wrapping point layers so that only they are rasterized
#'
#' @param raster NULL to rasterize above `threshold` cells, or TRUE / FALSE
#'
#' @param n_cells Number of cells plotted
#'
#' @param threshold Number of cells above which `raster = NULL` rasterizes
#'
#' @param dpi Resolution of the rasterized points
#'
#' @return A function taking and returning a layer
#'
#' @keywords internal
points_layer_fun <- function(raster,
                             n_cells,
                             threshold,
                             dpi) {
    do_raster <- if (is.null(raster)) n_cells > threshold else raster
    if (do_raster) check_installed("ggrastr", reason = "to rasterize points")
    function(layer) {
        if (do_raster) ggrastr::rasterise(layer, dpi = dpi) else layer
    }
}



#' Colour of the gap inside outlines and of label halos
#'
#' @param bg Background colour
#'
#' @keywords internal
halo_colour <- function(bg) {
    if (is.na(bg) || identical(bg, "transparent")) "white" else bg
}


#' Outline layers drawn beneath all cells
#'
#' A coloured layer then a background-coloured layer slightly larger than the
#' points, drawn beneath all cells, give the cell clouds a crisp border.
#'
#' @param df Plotted cells
#'
#' @param pt.size Point size
#'
#' @param outline.size Thickness of the outline
#'
#' @param outline.color A single colour, or one colour per row of `df`
#'
#' @param halo.color Colour of the gap between the outline and the cells
#'
#' @param layer_fun Function wrapping point layers, from
#'                  \code{\link{points_layer_fun}}
#'
#' @return A list of two layers
#'
#' @keywords internal
outline_layers <- function(df,
                           pt.size,
                           outline.size,
                           outline.color,
                           halo.color,
                           layer_fun) {
    outline <- if (length(outline.color) == 1) {
        geom_point(data = df, size = pt.size + outline.size, colour = outline.color, stroke = 0)
    } else {
        df$outline <- outline.color
        geom_point(aes(colour = I(outline)), data = df, size = pt.size + outline.size, stroke = 0)
    }
    list(
        layer_fun(outline),
        layer_fun(geom_point(data = df, size = pt.size + outline.size / 2, colour = halo.color, stroke = 0))
    )
}


#' Outlines around the cell clouds of each group, in data units
#'
#' Cells are converted to inches on the panel, outlined by
#' \code{\link{cloud_outline}} per group (and per panel, when split), and
#' converted back, so that the outlines keep the same distance from the cells
#' and the same smoothness on the page whatever the range of the embedding.
#'
#' @param df Cells with x, y, mark and, when `by_panel`, panel columns
#'
#' @param by_panel Whether to outline groups within each panel
#'
#' @param xlim,ylim Axis limits of the panel
#'
#' @param panel.width,panel.height Size of a panel in inches
#'
#' @inheritParams cloud_outline
#'
#' @return A data frame with x, y, mark, ring (one per closed outline) and,
#'         when `by_panel`, panel columns
#'
#' @keywords internal
mark_paths <- function(df,
                       by_panel,
                       xlim,
                       ylim,
                       panel.width,
                       panel.height,
                       expand,
                       smooth,
                       coverage) {
    df <- df[!is.na(df$mark), , drop = FALSE]
    x_per_in <- diff(xlim) / panel.width
    y_per_in <- diff(ylim) / panel.height
    sets <- if (by_panel) split(df, list(df$mark, df$panel), drop = TRUE) else split(df, df$mark, drop = TRUE)
    paths <- lapply(seq_along(sets), function(i) {
        d <- sets[[i]]
        out <- cloud_outline(
            (d$x - xlim[1]) / x_per_in, (d$y - ylim[1]) / y_per_in,
            expand = expand, smooth = smooth, coverage = coverage
        )
        if (nrow(out) == 0) return(NULL)
        res <- data.frame(
            x = xlim[1] + out$u * x_per_in, y = ylim[1] + out$v * y_per_in,
            mark = d$mark[1], ring = paste(i, out$ring, sep = ".")
        )
        if (by_panel) res$panel <- d$panel[1]
        res
    })
    do.call(rbind, paths)
}


#' Smooth outline around the cloud of a set of cells
#'
#' The cells are binned onto a grid and smoothed with a Gaussian kernel of sd
#' `smooth`. The region where this density is at least its value at the
#' `1 - coverage` quantile of the cells is the cloud: the sparsest cells, such
#' as strays far from the rest, fall outside it. Holes in the cloud are
#' filled. The cloud is then blurred and cut at the level that, along a
#' straight edge, lies exactly `expand` beyond it, which both pushes the
#' outline away from the cells and rounds its corners. Anything enclosed by
#' the outline is filled again so that only outer outlines are returned.
#'
#' @param u,v Coordinates of the cells, in inches
#'
#' @param expand Distance between the cloud and its outline, in inches
#'
#' @param smooth Smoothing of the cloud, in inches. Larger values give rounder
#'               outlines that bridge larger gaps between cells
#'
#' @param coverage Fraction of the cells, from the densest regions, taken to
#'                 form the cloud
#'
#' @param resolution Spacing of the grid, in inches
#'
#' @return A data frame with u, v and ring columns, one closed ring per
#'         separate part of the cloud
#'
#' @keywords internal
cloud_outline <- function(u,
                          v,
                          expand,
                          smooth,
                          coverage,
                          resolution = 0.01) {
    # Blurring with an sd of at least half the distance keeps the cut level
    # at no less than 2.3% of the blurred cloud, well clear of numerical noise
    blur <- max(smooth, expand / 2)
    edge <- stats::pnorm(-expand / blur)
    margin <- expand + 3 * blur + 3 * smooth
    gu <- seq(min(u) - margin, max(u) + margin + resolution, by = resolution)
    gv <- seq(min(v) - margin, max(v) + margin + resolution, by = resolution)
    kernel <- function(g, sd) stats::dnorm(outer(g, g, "-"), sd = sd) * resolution
    bu <- round((u - gu[1]) / resolution) + 1
    bv <- round((v - gv[1]) / resolution) + 1

    counts <- matrix(0, length(gu), length(gv))
    sums <- rowsum(rep(1, length(u)), bu + (bv - 1) * length(gu))
    counts[as.integer(rownames(sums))] <- sums
    density <- kernel(gu, smooth) %*% counts %*% kernel(gv, smooth)
    level <- stats::quantile(density[cbind(bu, bv)], 1 - coverage, names = FALSE)
    cloud <- fill_holes(density >= level)

    field <- kernel(gu, blur) %*% (cloud * 1) %*% kernel(gv, blur)
    field[fill_holes(field >= edge) & field < edge] <- 1
    iso <- isoband::isolines(gu, gv, t(field), edge)[[1]]
    data.frame(u = iso$x, v = iso$y, ring = iso$id)
}


#' Fill the holes of a logical mask
#'
#' @param mask A logical matrix
#'
#' @return The mask with every FALSE cell that cannot be reached from the
#'         border without crossing a TRUE cell set to TRUE
#'
#' @keywords internal
fill_holes <- function(mask) {
    open <- !mask
    n <- nrow(mask)
    m <- ncol(mask)
    reach <- matrix(FALSE, n, m)
    reach[c(1, n), ] <- open[c(1, n), ]
    reach[, c(1, m)] <- open[, c(1, m)]
    repeat {
        grown <- reach
        grown[-1, ] <- grown[-1, ] | reach[-n, ]
        grown[-n, ] <- grown[-n, ] | reach[-1, ]
        grown[, -1] <- grown[, -1] | reach[, -m]
        grown[, -m] <- grown[, -m] | reach[, -1]
        grown <- grown & open
        if (sum(grown) == sum(reach)) break
        reach <- grown
    }
    mask | (open & !reach)
}


#' Dashed outlines around the cell clouds of each group
#'
#' @param paths Output of \code{\link{mark_paths}}, with a colour column
#'
#' @param linetype,linewidth Line type and width of the outlines
#'
#' @return A polygon layer, drawn in vector form
#'
#' @keywords internal
mark_layer <- function(paths,
                       linetype,
                       linewidth) {
    geom_polygon(
        aes(x = x, y = y, group = ring, colour = I(colour)), data = paths,
        fill = NA, linetype = linetype, linewidth = linewidth, inherit.aes = FALSE,
        show.legend = FALSE
    )
}


#' Fixed axis limits without expansion
#'
#' @param rd Output of \code{\link{reduction_data}}
#'
#' @return A list of x and y scales
#'
#' @keywords internal
extent_scales <- function(rd) {
    list(
        scale_x_continuous(limits = rd$xlim, expand = expansion(0), oob = scales::oob_keep),
        scale_y_continuous(limits = rd$ylim, expand = expansion(0), oob = scales::oob_keep)
    )
}


#' Theme of the reduction plots: the embedding look with fixed panel sizes
#'
#' The embedding layout of \code{\link{theme_bioToolkit}}, fixed panel sizes
#' from \code{\link{panel_size_theme}}, and the legend and strip settings of
#' the reduction plots. Text sizes given as NULL keep the theme's.
#'
#' @param panel.width,panel.height Size of a single panel in inches
#'
#' @param ncol,nrow Number of panel columns and rows
#'
#' @param panel.spacing Gap between panels in inches
#'
#' @param legend.position Position of the legend
#'
#' @param legend.text.size Font size of the legend text
#'
#' @param strip.text.size,strip.text.face Font size and face of panel labels
#'
#' @param bg Background colour
#'
#' @param arrow.scale The `arrow.scale` of the dimension arrows, default is 1
#'
#' @return A ggplot2 theme
#'
#' @keywords internal
embedding_theme <- function(panel.width,
                            panel.height,
                            ncol,
                            nrow,
                            panel.spacing,
                            legend.position,
                            legend.text.size,
                            strip.text.size,
                            strip.text.face,
                            bg,
                            arrow.scale = 1) {
    # Room for the arrow labels below and left of the panels: the theme's
    # 0.2 cm plus 0.4 cm scaled as the arrows are (arrow_scale_factor()),
    # i.e. the theme's 0.6 cm at a 3.5-inch panel. The same in every plot of
    # a figure, with or without arrows, so that they line up
    k <- arrow_scale_factor(panel.width, panel.height, arrow.scale)
    arrow_room <- 0.2 + 0.4 * k
    theme_bioToolkit("embedding", bg = bg) +
        panel_size_theme(panel.width, panel.height, ncol = ncol, nrow = nrow,
                         panel.spacing = panel.spacing) +
        theme(
            legend.position = legend.position,
            # Properties left NULL inherit the theme's
            legend.text = element_text(size = legend.text.size),
            strip.text = element_text(size = strip.text.size, face = strip.text.face),
            plot.margin = margin(0.2, 0.2, arrow_room, arrow_room, "cm")
        )
}


#' Labels at the median position of each group
#'
#' @param df Cells with x, y, group and, when `by_panel`, panel columns
#'
#' @param by_panel Whether to label groups within each panel
#'
#' @param label.size Size of the labels
#'
#' @param label.color Colour of the labels
#'
#' @param halo.color Colour of the halo around the labels
#'
#' @param seed Seed for label placement
#'
#' @return A ggrepel text layer
#'
#' @keywords internal
group_label_layer <- function(df,
                              by_panel,
                              label.size,
                              label.color,
                              halo.color,
                              seed) {
    check_installed("ggrepel", reason = "to label groups")
    by <- if (by_panel) cbind(x, y) ~ group + panel else cbind(x, y) ~ group
    centers <- stats::aggregate(by, data = df, FUN = stats::median)
    ggrepel::geom_text_repel(
        data = centers,
        aes(label = group),
        size = label.size,
        fontface = "bold",
        colour = label.color,
        bg.colour = halo.color,
        bg.r = 0.15,
        seed = seed,
        min.segment.length = Inf
    )
}



#' Add dimension arrows to one or all panels
#'
#' @param p A ggplot object
#'
#' @param reduction Name of the reduction
#'
#' @param labels The two arrow labels
#'
#' @param arrow_args Further arguments to \code{\link{dimension_arrow_grob}}
#'
#' @param panel Levels of the facetting variable `panel` to draw the arrows
#'              in, or NULL to draw them in every panel
#'
#' @return a ggplot object
#'
#' @keywords internal
add_reduction_arrows <- function(p,
                                 reduction,
                                 labels,
                                 arrow_args,
                                 panel = NULL) {
    # Arrow labels take the theme's large text size unless given; the arrow
    # geometry keeps the defaults of dimension_arrow_grob()
    if (is.null(arrow_args$text_font_size)) arrow_args$text_font_size <- theme_sizes()$large
    # y_coord_adjust defaults to x_coord_adjust in dimension_arrow_grob()
    grob <- do.call(dimension_arrow_grob, c(list(reduction = reduction, dims = labels), arrow_args))
    if (is.null(panel)) {
        apply_arrow_grob(p, grob)
    } else if (length(panel) == 0) {
        # No panel gets arrows (`arrow.panels = "none"`); the axes are still
        # replaced by nothing, as in the panels without arrows
        apply_arrow_grob(p, grid::nullGrob())
    } else {
        lv <- levels(p$data$panel)
        apply_arrow_grob(p, grob, data = data.frame(panel = factor(panel, levels = lv, exclude = NULL)))
    }
}



