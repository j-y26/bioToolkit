# Some useful functions for enhancing the visualization of single cell data

#' @title Add dimension arrows
#'
#' @description Adds two short arrows with the dimension names (e.g. "UMAP-1",
#'              "UMAP-2") to the bottom left corner of an embedding plot and
#'              removes its axes, a compact replacement for full axes. Works
#'              on any ggplot, e.g. from \code{Seurat::DimPlot()}; the
#'              package's own embedding plots draw these arrows themselves.
#'
#'              The arrows keep a fixed proportion to the panel they are
#'              drawn in: the sizes below are those of a panel 3.5 inches on
#'              its shorter side, and everything (lengths, line width, label
#'              size) is scaled with the panel's actual size when the plot is
#'              drawn. A small panel gets small arrows that stay clear of the
#'              cells, a large one arrows that are not lost in it.
#'              `arrow.scale` makes them bigger or smaller at any panel size.
#'
#' @param plot A ggplot object
#'
#' @param reduction The reduction method used to generate the plot, default is
#'                 "umap", allowed values are "umap", "tsne", "pca", "spatial"
#'
#' @param dims A vector of dimension names, overwrites the default dimension
#'             names assumed from reduction method
#'
#' @param arrow.scale Size of the arrows relative to the panel, default is
#'                    1: arrows 2 cm long on a panel of 3.5 inches (about
#'                    23\% of its shorter side), scaled with the panel. Use
#'                    e.g. 0.8 for smaller or 1.25 for larger arrows
#'
#' @param arrow.length The length of the arrows in cm on a 3.5-inch panel,
#'                     default is 2
#'
#' @param x.coord.adjust Horizontal adjustment of the arrow origin in cm on a
#'                       3.5-inch panel, default is 0.25
#'
#' @param y.coord.adjust Vertical adjustment of the arrow origin, default is
#'                       the same as `x.coord.adjust`
#'
#' @param arrow.size The length of the arrow head in cm on a 3.5-inch panel,
#'                   default is 0.25
#'
#' @param arrow.width The width of the arrow head in cm on a 3.5-inch panel,
#'                    default is 0.2
#'
#' @param arrow.line.width The width of the arrow line on a 3.5-inch panel
#'                         (in R's `lwd` units), default is 3
#'
#' @param arrow.color The color of the arrow and label, default is "black"
#'
#' @param text.font.size The font size of the labels in points on a 3.5-inch
#'                       panel, default is 12
#'
#' @param text.dst The distance of the labels from the arrows in cm on a
#'                 3.5-inch panel, default is 0.3
#'
#' @param ... Deprecated snake_case names of the arguments above (e.g.
#'            `arrow_length`), accepted with a warning; renamed to dot.case
#'            in 0.10.0
#'
#' @return a ggplot object with the dimension arrows added
#'
#' @import ggplot2
#' @import grid
#'
#' @export
#'
add_dimension_arrows <- function(plot,
                                 reduction = "umap",
                                 dims = NULL,
                                 arrow.scale = 1,
                                 arrow.length = 2,
                                 x.coord.adjust = 0.25,
                                 y.coord.adjust = x.coord.adjust,
                                 arrow.size = 0.25,
                                 arrow.width = 0.2,
                                 arrow.line.width = 3,
                                 arrow.color = "black",
                                 text.font.size = 12,
                                 text.dst = 0.3,
                                 ...) {
  # Old snake_case names (before 0.10.0) are mapped to the new ones, with a
  # deprecation warning; y.coord.adjust's default is evaluated only after this
  list2env(renamed_args(list(...), ARROW_ARG_RENAMES, "add_dimension_arrows",
                        names(match.call())), environment())

  # Check if the plot is a ggplot object
  if (!is_ggplot(plot)) {
    stop("The plot must be a ggplot object")
  }

  # The internal grob builder keeps its own (snake_case) argument names
  arrow_grob <- dimension_arrow_grob(
    reduction = reduction,
    dims = dims,
    arrow_scale = arrow.scale,
    arrow_length = arrow.length,
    x_coord_adjust = x.coord.adjust,
    y_coord_adjust = y.coord.adjust,
    arrow_size = arrow.size,
    arrow_width = arrow.width,
    arrow_line_width = arrow.line.width,
    arrow_color = arrow.color,
    text_font_size = text.font.size,
    text_dst = text.dst
  )

  apply_arrow_grob(plot, arrow_grob)
}


# Old (before 0.10.0) -> new argument names of add_dimension_arrows(), also
# accepted by format_sc_plots()
ARROW_ARG_RENAMES <- c(
  arrow_length = "arrow.length",
  x_coord_adjust = "x.coord.adjust",
  y_coord_adjust = "y.coord.adjust",
  arrow_size = "arrow.size",
  arrow_width = "arrow.width",
  arrow_line_width = "arrow.line.width",
  arrow_color = "arrow.color",
  text_font_size = "text.font.size",
  text_dst = "text.dst"
)


# Dimension arrows scaled with the panel (revision of 2026-10-15)
# ---------------------------------------------------------------
# The arrows used to have a fixed size in cm: right on a 3.5-4 inch panel,
# but on a small panel they were too big and their corner, at a fixed 0.25 cm
# from the panel edge, ran into cells near the left or bottom edge (the
# panel's 5% expansion is only 0.19 cm on a 1.5-inch panel). User decision
# (2026-10-15): the arrows keep the proportions they have on a 3.5-inch panel
# (ARROW_REF_PANEL) at every panel size. All of their geometry (lengths, line
# width, label size) is multiplied by
#     arrow.scale * (shorter side of the panel) / 3.5 in
# (arrow_scale_factor()). On a 3.5-inch panel the corner and arrow heads lie
# within 4% of the edges, inside the 5% expansion of the axis limits, so
# they stay clear of the cells at any size. The scaling is done when the
# grob is drawn (makeContent()), from the size of the panel it is drawn in,
# so it also works for plots without a fixed panel size (e.g. Seurat plots
# passed to format_sc_plots()) and follows the window when one is resized.

# Shorter side of the panel (inches) at which the arrows have the sizes given
# by the arguments of add_dimension_arrows()
ARROW_REF_PANEL <- 3.5


#' Factor by which the dimension arrows are scaled for a panel
#'
#' @param width,height Size of the panel in inches
#'
#' @param arrow_scale The `arrow.scale` given by the user
#'
#' @return `arrow_scale` times the panel's shorter side over
#'         `ARROW_REF_PANEL` (3.5 inches)
#'
#' @keywords internal
arrow_scale_factor <- function(width,
                               height,
                               arrow_scale = 1) {
  arrow_scale * min(width, height) / ARROW_REF_PANEL
}


#' Check the `arrow.scale` argument
#'
#' @param arrow_scale The value given
#'
#' @return Nothing; stops when the value is not a single positive number
#'
#' @keywords internal
check_arrow_scale <- function(arrow_scale) {
  if (!(is.numeric(arrow_scale) && length(arrow_scale) == 1 && !is.na(arrow_scale) && arrow_scale > 0)) {
    stop("`arrow.scale` must be a single positive number", call. = FALSE)
  }
}


#' Build the dimension arrow grob
#'
#' Returns an empty grob tree of class "bt_dimension_arrows" that holds the
#' arrow geometry; the arrows and labels themselves are built when it is
#' drawn (\code{makeContent.bt_dimension_arrows()}), scaled with the size of
#' the viewport (i.e. panel) it is drawn in, and measured from its bottom left
#' corner.
#'
#' @param reduction,dims As in \code{\link{add_dimension_arrows}}
#'
#' @param arrow_scale,arrow_length,x_coord_adjust,y_coord_adjust,arrow_size,arrow_width,arrow_line_width,arrow_color,text_font_size,text_dst
#'        The arrow geometry, as the dot.case arguments of
#'        \code{\link{add_dimension_arrows}} (e.g. `arrow_length` is
#'        `arrow.length`): sizes on a 3.5-inch panel
#'
#' @return A grob tree
#'
#' @keywords internal
dimension_arrow_grob <- function(reduction = "umap",
                                 dims = NULL,
                                 arrow_scale = 1,
                                 arrow_length = 2,
                                 x_coord_adjust = 0.25,
                                 y_coord_adjust = x_coord_adjust,
                                 arrow_size = 0.25,
                                 arrow_width = 0.2,
                                 arrow_line_width = 3,
                                 arrow_color = "black",
                                 text_font_size = 12,
                                 text_dst = 0.3) {
  check_arrow_scale(arrow_scale)
  # Get the dimension names if not provided
  if (is.null(dims)) {
    dims <- reduction_dim_labels(reduction)
  }

  # The geometry is kept as given (sizes on a 3.5-inch panel); tagged so that
  # the arrows of a plot can be found and replaced (remove_dimension_arrows())
  gTree(
    dims = dims,
    arrow_scale = arrow_scale,
    geometry = list(
      arrow_length = arrow_length,
      x_coord_adjust = x_coord_adjust,
      y_coord_adjust = y_coord_adjust,
      arrow_size = arrow_size,
      arrow_width = arrow_width,
      arrow_line_width = arrow_line_width,
      text_font_size = text_font_size,
      text_dst = text_dst
    ),
    arrow_color = arrow_color,
    cl = "bt_dimension_arrows"
  )
}


#' @importFrom grid makeContent
#' @export
makeContent.bt_dimension_arrows <- function(x) {
  # Size of the panel the arrows are drawn in, and the scale factor for it
  width <- convertWidth(unit(1, "npc"), "in", valueOnly = TRUE)
  height <- convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)
  k <- arrow_scale_factor(width, height, x$arrow_scale)
  # Every size (cm, lwd, points) is scaled by the same factor
  g <- lapply(x$geometry, function(v) v * k)
  setChildren(x, dimension_arrow_children(x$dims, g, x$arrow_color))
}


#' The arrows and labels of the dimension arrow grob
#'
#' @param dims The two labels
#'
#' @param g The geometry, already scaled for the panel: a list with
#'          `arrow_length`, `x_coord_adjust`, `y_coord_adjust`, `arrow_size`,
#'          `arrow_width` and `text_dst` in cm, `arrow_line_width` (lwd) and
#'          `text_font_size` (points)
#'
#' @param arrow_color Colour of the arrows and labels
#'
#' @return A gList
#'
#' @keywords internal
dimension_arrow_children <- function(dims,
                                     g,
                                     arrow_color) {
  cm <- function(v) unit(v, "cm")
  x0 <- g$x_coord_adjust
  y0 <- g$y_coord_adjust
  # Where the shafts end and the heads begin
  x_head <- x0 + g$arrow_length - g$arrow_size
  y_head <- y0 + g$arrow_length - g$arrow_size

  # Arrow shaft grobs (without arrowheads)
  arrow_x <- segmentsGrob(
    x0 = cm(x0),
    y0 = cm(y0),
    x1 = cm(x_head),
    y1 = cm(y0),
    gp = gpar(col = arrow_color, lwd = g$arrow_line_width)
  )

  arrow_y <- segmentsGrob(
    x0 = cm(x0),
    y0 = cm(y0),
    x1 = cm(x0),
    y1 = cm(y_head),
    gp = gpar(col = arrow_color, lwd = g$arrow_line_width)
  )

  # Arrowhead polygons (triangles) at end of arrows
  arrowhead_x <- polygonGrob(
    x = cm(c(x_head, x0 + g$arrow_length, x_head)),
    y = cm(c(y0 + g$arrow_width / 2, y0, y0 - g$arrow_width / 2)),
    gp = gpar(fill = arrow_color, col = arrow_color)
  )

  arrowhead_y <- polygonGrob(
    x = cm(c(x0 - g$arrow_width / 2, x0, x0 + g$arrow_width / 2)),
    y = cm(c(y_head, y0 + g$arrow_length, y_head)),
    gp = gpar(fill = arrow_color, col = arrow_color)
  )

  # Labels for the arrows, below the horizontal one and left of the vertical
  # one
  label_gp <- gpar(col = arrow_color, fontsize = g$text_font_size, fontface = "bold")
  label_x <- textGrob(
    label = dims[1],
    x = cm(x0 + g$arrow_length / 2),
    y = cm(y0 - g$text_dst),
    just = "center",
    gp = label_gp
  )

  label_y <- textGrob(
    label = dims[2],
    x = cm(x0 - g$text_dst),
    y = cm(y0 + g$arrow_length / 2),
    rot = 90,
    just = "center",
    gp = label_gp
  )

  gList(arrow_x, arrow_y, arrowhead_x, arrowhead_y, label_x, label_y)
}


#' Axis labels for the dimensions of a reduction
#'
#' @param reduction Name of the reduction, e.g. "umap"
#'
#' @param dims The two dimensions plotted, default is 1 and 2
#'
#' @return A character vector of two labels, e.g. "UMAP-1" and "UMAP-2"
#'
#' @keywords internal
reduction_dim_labels <- function(reduction,
                                 dims = c(1, 2)) {
  reduction <- tolower(reduction)
  prefix <- if (grepl("umap", reduction)) {
    "UMAP"
  } else if (grepl("tsne", reduction) | grepl("t-sne", reduction)) {
    "tSNE"
  } else if (grepl("pca", reduction)) {
    "PC"
  } else if (grepl("spatial", reduction)) {
    "Spatial"
  } else {
    reduction
  }
  paste0(prefix, "-", dims)
}


#' Draw an arrow grob in one or all panels of a plot
#'
#' `annotation_custom()` has no data of its own and is therefore drawn in every
#' panel. Passing `data` a single row that carries the facetting variable
#' restricts the grob to the matching panel, which is how the arrows are kept to
#' the bottom left panel of a split plot.
#'
#' @param plot A ggplot object
#'
#' @param grob The grob to draw, e.g. from \code{\link{dimension_arrow_grob}}
#'
#' @param data A one-row data frame naming the panel to draw in, or NULL
#'             (default) to draw in every panel
#'
#' @return a ggplot object
#'
#' @keywords internal
apply_arrow_grob <- function(plot,
                             grob,
                             data = NULL) {
  annotation <- if (is.null(data)) {
    annotation_custom(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
  } else {
    layer(
      geom = GeomCustomAnn, stat = StatIdentity, position = PositionIdentity,
      data = data, mapping = NULL, inherit.aes = FALSE, show.legend = FALSE,
      params = list(grob = grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
    )
  }

  # The arrows are drawn outside the panel, so clipping is switched off. The
  # coordinate system is replaced directly (adding a second one makes
  # ggplot2 replace it with a message), keeping the limits of a cartesian one
  coords <- plot$coordinates
  if (!identical(coords$clip, "off")) {
    plot$coordinates <- if (inherits(coords, "CoordCartesian") && !isTRUE(coords$default)) {
      coord_cartesian(xlim = coords$limits$x, ylim = coords$limits$y, expand = coords$expand, clip = "off")
    } else {
      coord_cartesian(clip = "off")
    }
  }
  plot +
    annotation +
    theme(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank()
    )
}

#' Remove the dimension arrows of a plot
#'
#' @param plot A ggplot
#'
#' @return The plot without the layers drawing \code{dimension_arrow_grob()}
#'         arrows
#'
#' @keywords internal
remove_dimension_arrows <- function(plot) {
  is_arrow <- vapply(plot$layers, function(l) {
    inherits(l$geom_params$grob, "bt_dimension_arrows")
  }, logical(1))
  plot$layers <- plot$layers[!is_arrow]
  plot
}


#' @title Restyle a single ggplot made by Seurat
#'
#' @description Removes the axes, legend or grid of one ggplot (e.g. from
#'              \code{Seurat::DimPlot()}), restores a hidden legend, or sets
#'              the legend position and title. Used by
#'              \code{\link{format_sc_plots}} on each plot of a list; only
#'              valid for individual ggplot objects.
#'
#' @param plot A single ggplot object
#'
#' @param remove.axes Whether to remove the axes, default is FALSE
#'
#' @param remove.legend Whether to remove the legend, default is FALSE
#'
#' @param remove.grid Whether to remove the grid, default is FALSE
#'
#' @param restore.legend Whether to restore the legend, default is FALSE
#'
#' @param legend.position The position of the legend, default is NULL (meaning
#'                        no change)
#'
#' @param legend.title The title of the legend, default is NULL
#'
#' @param legend.title.size The size of the legend title, default is 12
#'
#' @param legend.title.bold Whether to bold the legend title, default is TRUE
#'
#' @param ... Deprecated snake_case names of the arguments above (e.g.
#'            `remove_axes`), accepted with a warning; renamed to dot.case
#'            in 0.10.0
#'
#' @import ggplot2
#' @export
#'
#' @return a ggplot object with the axes, legend, and grid manipulated
#'
manipulate_sc_plot <- function(plot,
                               remove.axes = FALSE,
                               remove.legend = FALSE,
                               remove.grid = FALSE,
                               restore.legend = FALSE,
                               legend.position = NULL,
                               legend.title = NULL,
                               legend.title.size = 12,
                               legend.title.bold = TRUE,
                               ...) {
  # Old snake_case names (before 0.10.0) are mapped to the new ones
  list2env(renamed_args(list(...), PLOT_ARG_RENAMES, "manipulate_sc_plot",
                        names(match.call())), environment())

  # Check if the plot is a ggplot object
  if (!is_ggplot(plot)) {
    stop("The plot must be a ggplot object")
  }

  # The Seurat helpers (NoAxes() etc.) are plain theme changes; Seurat is in
  # Suggests
  if (remove.axes || remove.legend || remove.grid || restore.legend) {
    check_installed("Seurat", reason = "to restyle Seurat plots")
  }
  # Remove axes
  if (remove.axes) {
    plot <- plot + Seurat::NoAxes()
  }
  # Remove legend
  if (remove.legend) {
    plot <- plot + Seurat::NoLegend()
  }
  # Remove grid
  if (remove.grid) {
    plot <- plot + Seurat::NoGrid()
  }
  # Restore legend
  if (restore.legend) {
    plot <- plot + Seurat::RestoreLegend()
  }
  # Set legend position
  if (!is.null(legend.position)) {
    plot <- plot + theme(legend.position = legend.position)
  }
  # Set legend title
  if (!is.null(legend.title)) {
    plot <- plot + labs(color = legend.title)
    face <- if (legend.title.bold) "bold" else "plain"
    plot <- plot + theme(legend.title = element_text(face = face, size = legend.title.size))
  }

  return(plot)
}


# Old (before 0.10.0) -> new argument names of manipulate_sc_plot(), also
# accepted by format_sc_plots()
PLOT_ARG_RENAMES <- c(
  remove_axes = "remove.axes",
  remove_legend = "remove.legend",
  remove_grid = "remove.grid",
  restore_legend = "restore.legend",
  legend_position = "legend.position",
  legend_title = "legend.title",
  legend_title_size = "legend.title.size",
  legend_title_bold = "legend.title.bold"
)


#' @title Format a set of embedding plots (from Seurat, scCustomize, SCpubr
#'        etc.)
#'
#' @description Formats a list of ggplot objects or a patchwork object
#'              consistently and arranges them in a grid. By default the axes
#'              are removed and replaced by dimension arrows in the bottom left
#'              corner of the bottom left plot; `arrow.panels` puts them on
#'              every plot or on one edge of the grid instead. Arrows a plot
#'              already has (e.g. from \code{\link{DimPlot_bioToolkit}}) are
#'              replaced, so that each plot has at most one set.
#'
#' @param plots A list of ggplot objects or a patchwork object
#'
#' @param ncol Number of columns in the plot grid, default is NULL. Prioritized
#'             if nrow is also set.
#'
#' @param nrow Number of rows in the plot grid, default is NULL
#'
#' @param by.row Whether to fill the plot grid by row, default is TRUE,
#'               consistent with patchwork's default behavior. If FALSE, the
#'               plots will be filled by column.
#'
#' @param arrow.panels Plots that get the dimension arrows: "bottom.left"
#'                     (default; the first plot of the bottom row), "all",
#'                     "left" / "right" (the first / last plot of each row),
#'                     "top" / "bottom" (the top / lowest plot of each column)
#'                     or "none". Edges are those of the plots drawn: in a grid
#'                     whose last row is incomplete, a plot with nothing below
#'                     it is on the bottom edge. The grid is given by `ncol` /
#'                     `nrow`, the patchwork's layout, or patchwork's default.
#'                     Unless "none", the axes of all plots are removed.
#'                     (Replaces `add.arrow` and `bottom.left.only`, accepted
#'                     with a deprecation warning)
#'
#' @param remove.axes Whether to remove the axes, default is TRUE
#'
#' @param reduction The reduction method used to generate the plots, default
#'                  is "umap", allowed values are "umap", "tsne", "pca",
#'                  "spatial". Used only when arrows are drawn
#'
#' @param dims A vector of dimension names, overwrites the default dimension
#'             names. Used only when arrows are drawn
#'
#' @param ... Deprecated arguments, accepted with a warning: `add.arrow` and
#'            `bottom.left.only` (replaced by `arrow.panels`), and the
#'            snake_case names renamed to dot.case in 0.10.0 (e.g. `byrow`,
#'            `legend_position`, `arrow_length`)
#'
#' @inheritParams manipulate_sc_plot
#' @inheritParams add_dimension_arrows
#'
#' @import ggplot2
#' @import patchwork
#'
#' @export
#'
#' @return a patchwork object with the formatted plots
#'
format_sc_plots <- function(plots,
                            ncol = NULL,
                            nrow = NULL,
                            by.row = TRUE,
                            arrow.panels = c("bottom.left", "all", "left", "right", "top", "bottom", "none"),
                            remove.axes = TRUE,
                            remove.legend = FALSE,
                            remove.grid = FALSE,
                            restore.legend = FALSE,
                            legend.position = NULL,
                            legend.title = NULL,
                            legend.title.size = 12,
                            legend.title.bold = TRUE,
                            reduction = "umap",
                            dims = NULL,
                            arrow.scale = 1,
                            arrow.length = 2,
                            x.coord.adjust = 0.25,
                            y.coord.adjust = x.coord.adjust,
                            arrow.size = 0.25,
                            arrow.width = 0.2,
                            arrow.line.width = 3,
                            arrow.color = "black",
                            text.font.size = 12,
                            text.dst = 0.3,
                            ...) {
  arrow.panels <- match.arg(arrow.panels)
  dots <- list(...)
  # `add.arrow` and `bottom.left.only` (before 1.0.0; also as add_arrow /
  # bottom_left_only) are translated into `arrow.panels`
  old_arrow <- intersect(names(dots), c("add.arrow", "bottom.left.only", "add_arrow", "bottom_left_only"))
  if (length(old_arrow)) {
    if ("arrow.panels" %in% names(match.call())) {
      stop("In `format_sc_plots()`, give either `arrow.panels` or the deprecated `",
           old_arrow[1], "`, not both", call. = FALSE)
    }
    warning(warningCondition(
      paste0("In `format_sc_plots()`, ", paste0("`", old_arrow, "`", collapse = ", "),
             " deprecated; use `arrow.panels` instead."),
      class = "deprecatedWarning"
    ))
    add <- dots[["add.arrow"]] %||% dots[["add_arrow"]] %||% TRUE
    bl <- dots[["bottom.left.only"]] %||% dots[["bottom_left_only"]] %||% FALSE
    arrow.panels <- if (!isTRUE(add)) "none" else if (isTRUE(bl)) "bottom.left" else "all"
    dots <- dots[setdiff(names(dots), old_arrow)]
  }
  # Old snake_case names (before 0.10.0) are mapped to the new ones
  renames <- c(byrow = "by.row", PLOT_ARG_RENAMES, ARROW_ARG_RENAMES)
  list2env(renamed_args(dots, renames, "format_sc_plots", names(match.call())), environment())

  # Check if the plots is a list of ggplot objects or a patchwork object
  if (!is.list(plots) && !inherits(plots, "patchwork")) {
    stop("The plots must be a list of ggplot objects or a patchwork object")
  }

  # Check if the plots is a list of ggplot objects
  if (is.list(plots)) {
    for (i in seq_along(plots)) {
      if (!is_ggplot(plots[[i]])) {
        stop(paste0("The plot ", i, " must be a ggplot object"))
      }
    }
  }

  # Convert the patchwork object to a list of ggplot objects
  if (inherits(plots, "patchwork")) {
    plot_list <- as.list(plots)
  } else {
    plot_list <- plots
  }

  # Apply individual plot formatting
  plot_list <- lapply(plot_list, function(plot) {
    manipulate_sc_plot(
      plot = plot,
      remove.axes = remove.axes,
      remove.legend = remove.legend,
      remove.grid = remove.grid,
      restore.legend = restore.legend,
      legend.position = legend.position,
      legend.title = legend.title,
      legend.title.size = legend.title.size,
      legend.title.bold = legend.title.bold
    )
  })

  # Arrows with the same geometry on every plot that gets them
  add_arrows <- function(plot) {
    add_dimension_arrows(
      plot = plot,
      reduction = reduction,
      dims = dims,
      arrow.scale = arrow.scale,
      arrow.length = arrow.length,
      x.coord.adjust = x.coord.adjust,
      y.coord.adjust = y.coord.adjust,
      arrow.size = arrow.size,
      arrow.width = arrow.width,
      arrow.line.width = arrow.line.width,
      arrow.color = arrow.color,
      text.font.size = text.font.size,
      text.dst = text.dst
    )
  }

  # Layout of the grid: from ncol / nrow, else the patchwork's, else
  # patchwork's default (grDevices::n2mfrow(), as patchwork::wrap_dims())
  n_plots <- length(plot_list)
  if (is.null(ncol) && is.null(nrow) && inherits(plots, "patchwork")) {
    ncol <- plots$patches$layout$ncol
    nrow <- plots$patches$layout$nrow
    if (!is.null(plots$patches$layout$byrow)) by.row <- plots$patches$layout$byrow
  }
  if (!is.null(ncol)) {
    # ncol is used when both are given
    nrow <- ceiling(n_plots / ncol)
  } else if (!is.null(nrow)) {
    ncol <- ceiling(n_plots / nrow)
  } else {
    dims <- grDevices::n2mfrow(n_plots)
    nrow <- dims[2]
    ncol <- dims[1]
  }
  if (inherits(plots, "patchwork")) {
    plots$patches$layout$ncol <- ncol
    plots$patches$layout$nrow <- nrow
    plots$patches$layout$byrow <- by.row
  }

  # Arrows: replace any the plots already have, on the plots chosen by
  # `arrow.panels`; the axes of the others are removed too
  if (arrow.panels != "none") {
    idx <- arrow_panel_index(n_plots, ncol, nrow, by.row, arrow.panels)
    plot_list <- lapply(seq_len(n_plots), function(i) {
      plot <- remove_dimension_arrows(plot_list[[i]])
      if (i %in% idx) add_arrows(plot) else apply_arrow_grob(plot, grid::nullGrob())
    })
  }

  # Combine the plots into a patchwork object
  if (inherits(plots, "patchwork")) {
    combined_plot <- wrap_plots(plotlist = plot_list,
                                ncol = if (is.null(ncol)) NULL else ncol,
                                nrow = if (is.null(nrow)) NULL else nrow,
                                byrow = by.row,
                                widths = plots$patches$layout$widths,
                                heights = plots$patches$layout$heights,
                                guides = plots$patches$layout$guides,
                                tag_level = plots$patches$layout$tag_level,
                                design = plots$patches$layout$design,
                                axes = plots$patches$layout$axes,
                                axis_titles = plots$patches$layout$axis_titles)
    combined_plot <- combined_plot + plots$patches$annotation

  } else {
    combined_plot <- patchwork::wrap_plots(plotlist = plot_list, ncol = ncol, byrow = by.row)
  }

  return(combined_plot)
}
