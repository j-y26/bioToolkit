# Some useful functions for enhancing the visualization of single cell data

#' @title Add dimension arrows
#'
#' @description This function adds dimension arrows to a plot at the bottom
#'              left corner of the plot. The arrows represent the direction of
#'              the dimensions in the plot. Assumes that the plot is a ggplot
#'              object and does not have an existing axis.
#'
#' @param plot A ggplot object
#'
#' @param reduction The reduction method used to generate the plot, default is
#'                 "umap", allowed values are "umap", "tsne", "pca", "spatial"
#'
#' @param dims A vector of dimension names, overwrites the default dimension
#'             names assumed from reduction method
#'
#' @param arrow_length The size of the arrow, default is 2 (cm)
#'
#' @param x_coord_adjust Horizontal adjustment of the arrow origin, default is 0.25 (cm)
#'
#' @param y_coord_adjust Vertical adjustment of the arrow origin, default is the same as x_coord_adjust
#'
#' @param arrow_size The size of the arrow head, default is 0.25 cm
#'
#' @param arrow_width The width of the arrow, default is 0.2 cm
#'
#' @param arrow_line_width The width of the arrow line, default is 3
#'
#' @param arrow_color The color of the arrow and label, default is "black"
#'
#' @param text_font_size The font size of the text, default is 12
#'
#' @param text_dst The distance of the text from the arrow in cm, default is 0.3
#'
#'
#' @return a ggplot object with the dimension arrows added
#'
#'
#' @import ggplot2
#' @import grid
#'
#' @export
#'
add_dimension_arrows <- function(plot,
                                 reduction = "umap",
                                 dims = NULL,
                                 arrow_length = 2,
                                 x_coord_adjust = 0.25,
                                 y_coord_adjust = x_coord_adjust,
                                 arrow_size = 0.25,
                                 arrow_width = 0.2,
                                 arrow_line_width = 3,
                                 arrow_color = "black",
                                 text_font_size = 12,
                                 text_dst = 0.3) {
  # Check if the plot is a ggplot object
  if (!is_ggplot(plot)) {
    stop("The plot must be a ggplot object")
  }

  arrow_grob <- dimension_arrow_grob(
    reduction = reduction,
    dims = dims,
    arrow_length = arrow_length,
    x_coord_adjust = x_coord_adjust,
    y_coord_adjust = y_coord_adjust,
    arrow_size = arrow_size,
    arrow_width = arrow_width,
    arrow_line_width = arrow_line_width,
    arrow_color = arrow_color,
    text_font_size = text_font_size,
    text_dst = text_dst
  )

  apply_arrow_grob(plot, arrow_grob)
}


#' Build the dimension arrow grob
#'
#' Draws the two arrows and their labels in absolute units, measured from the
#' bottom left corner of whichever viewport (i.e. panel) the grob is drawn in.
#'
#' @inheritParams add_dimension_arrows
#'
#' @return A grob tree
#'
#' @keywords internal
dimension_arrow_grob <- function(reduction = "umap",
                                 dims = NULL,
                                 arrow_length = 2,
                                 x_coord_adjust = 0.25,
                                 y_coord_adjust = x_coord_adjust,
                                 arrow_size = 0.25,
                                 arrow_width = 0.2,
                                 arrow_line_width = 3,
                                 arrow_color = "black",
                                 text_font_size = 12,
                                 text_dst = 0.3) {
  # Get the dimension names if not provided
  if (is.null(dims)) {
    dims <- reduction_dim_labels(reduction)
  }

  # Arrow shaft grobs (without arrowheads)
  arrow_x <- segmentsGrob(
    x0 = unit(x_coord_adjust, "cm"), y0 = unit(y_coord_adjust, "cm"),
    x1 = unit(x_coord_adjust + arrow_length - arrow_size, "cm"), y1 = unit(y_coord_adjust, "cm"),
    gp = gpar(col = arrow_color, lwd = arrow_line_width)
  )

  arrow_y <- segmentsGrob(
    x0 = unit(x_coord_adjust, "cm"), y0 = unit(y_coord_adjust, "cm"),
    x1 = unit(x_coord_adjust, "cm"), y1 = unit(y_coord_adjust + arrow_length - arrow_size, "cm"),
    gp = gpar(col = arrow_color, lwd = arrow_line_width)
  )

  # Arrowhead polygons (triangles) at end of arrows
  arrowhead_x <- polygonGrob(
    x = unit.c(
      unit(x_coord_adjust + arrow_length - arrow_size, "cm"),
      unit(x_coord_adjust + arrow_length, "cm"),
      unit(x_coord_adjust + arrow_length - arrow_size, "cm")
    ),
    y = unit.c(
      unit(y_coord_adjust + arrow_width / 2, "cm"),
      unit(y_coord_adjust, "cm"),
      unit(y_coord_adjust - arrow_width / 2, "cm")
    ),
    gp = gpar(fill = arrow_color, col = arrow_color)
  )

  arrowhead_y <- polygonGrob(
    x = unit.c(
      unit(x_coord_adjust - arrow_width / 2, "cm"),
      unit(x_coord_adjust, "cm"),
      unit(x_coord_adjust + arrow_width / 2, "cm")
    ),
    y = unit.c(
      unit(y_coord_adjust + arrow_length - arrow_size, "cm"),
      unit(y_coord_adjust + arrow_length, "cm"),
      unit(y_coord_adjust + arrow_length - arrow_size, "cm")
    ),
    gp = gpar(fill = arrow_color, col = arrow_color)
  )

  # Labels for the arrows
  label_x <- textGrob(
    label = dims[1],
    x = unit(x_coord_adjust + arrow_length / 2, "cm"),
    y = unit(y_coord_adjust - text_dst, "cm"),
    just = "center",
    gp = gpar(col = arrow_color, fontsize = text_font_size, fontface = "bold")
  )

  label_y <- textGrob(
    label = dims[2],
    x = unit(x_coord_adjust - text_dst, "cm"),
    y = unit(y_coord_adjust + arrow_length / 2, "cm"),
    rot = 90,
    just = "center",
    gp = gpar(col = arrow_color, fontsize = text_font_size, fontface = "bold")
  )

  # Combine the grobs into a single grob
  grobTree(arrow_x, arrow_y, arrowhead_x, arrowhead_y, label_x, label_y)
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

  plot +
    annotation +
    coord_cartesian(clip = "off") +
    theme(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank()
    )
}

#' @title Helper to manipulate ggplots using Seurat's patchwork functions
#'
#' @description Only valid for individual ggplot objects
#'
#' @param plot A single ggplot object
#'
#' @param remove_axes Whether to remove the axes, default is FALSE
#'
#' @param remove_legend Whether to remove the legend, default is FALSE
#'
#' @param remove_grid Whether to remove the grid, default is FALSE
#'
#' @param restore_legend Whether to restore the legend, default is FALSE
#'
#' @param legend_position The position of the legend, default is NULL (meaning no change)
#'
#' @param legend_title The title of the legend, default is NULL
#'
#' @param legend_title_size The size of the legend title, default is 12
#'
#' @param legend_title_bold Whether to bold the legend title, default is FALSE
#'
#' @import ggplot2
#' @export
#'
#' @return a ggplot object with the axes, legend, and grid manipulated
#'
manipulate_sc_plot <- function(plot,
                     remove_axes = FALSE,
                     remove_legend = FALSE,
                     remove_grid = FALSE,
                     restore_legend = FALSE,
                     legend_position = NULL,
                     legend_title = NULL,
                     legend_title_size = 12,
                     legend_title_bold = FALSE) {
  # Check if the plot is a ggplot object
  if (!is_ggplot(plot)) {
    stop("The plot must be a ggplot object")
  }

  # Remove axes
  if (remove_axes) {
    if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required.")
    plot <- plot + Seurat::NoAxes()
  }
  # Remove legend
  if (remove_legend) {
    if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required.")
    plot <- plot + Seurat::NoLegend()
  }
  # Remove grid
  if (remove_grid) {
    if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required.")
    plot <- plot + Seurat::NoGrid()
  }
  # Restore legend
  if (restore_legend) {
    if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required.")
    plot <- plot + Seurat::RestoreLegend()
  }
  # Set legend position
  if (!is.null(legend_position)) {
    plot <- plot + theme(legend.position = legend_position)
  }
  # Set legend title
  if (!is.null(legend_title)) {
    plot <- plot + labs(color = legend_title)

    if (legend_title_bold) {
      plot <- plot + theme(legend.title = element_text(face = "bold", size = legend_title_size))
    } else {
      plot <- plot + theme(legend.title = element_text(size = legend_title_size))
    }
  }

  return(plot)
}


#' @title Format Dimension Plots (from Seurat, scCustomize, SCpubr etc.)
#'
#' @description This function operates on a list of
#'              ggplot objects or a patchwork object and formats them to have a
#'              consistent theme. The function performs group operations to
#'              format the plots. The default behavior is to remove the axes and
#'              add an arrow to the bottom left corner of each individual plot.
#'
#' @param plots A list of ggplot objects or a patchwork object
#'
#' @param ncol Number of columns in the plot grid, default is NULL. Prioritized if
#'             nrow is also set.
#'
#' @param nrow Number of rows in the plot grid, default is NULL
#'
#' @param byrow Whether to fill the plot grid by row, default is TRUE, consistent
#'              with patchwork's default behavior. If FALSE, the plots will be
#'              filled by column.
#'
#' @param add_arrow Whether to add an arrow to replace the axes, default is TRUE.
#'                  This will override the remove_axes argument and remove the axes
#'                  from all plots
#'
#' @param bottom_left_only Whether to add the arrow only to the plot at the bottom
#'                         left of the entire grid, default is FALSE. Only valid
#'                         if the arrangement of the plots are defined, either by
#'                         setting ncol or nrow, or that patchwork already defines
#'                         the layout. Only used when add_arrow is TRUE
#'
#' @param remove_axes Whether to remove the axes, default is TRUE
#'
#' @param remove_legend Whether to remove the legend, default is FALSE
#'
#' @param remove_grid Whether to remove the grid, default is FALSE
#'
#' @param restore_legend Whether to restore the legend, default is FALSE
#'
#' @param legend_position The position of the legend, default is NULL (meaning no change)
#'
#' @param legend_title The title of the legend, default is NULL
#'
#' @param legend_title_size The size of the legend title, default is 12
#'
#' @param legend_title_bold Whether to bold the legend title, default is FALSE
#'
#' @param reduction The reduction method used to generate the plots, default is
#'                 "umap", allowed values are "umap", "tsne", "pca", "spatial".
#'                 Used only if add_arrow is TRUE
#'
#' @param dims A vector of dimension names, overwrites the default dimension.
#'             Used only if add_arrow is TRUE
#'
#' @param arrow_length The size of the arrow, default is 2 (cm)
#'                     Used only if add_arrow is TRUE
#'
#' @param x_coord_adjust Horizontal adjustment of the arrow origin, default is 0.25 (cm)
#'                       Used only if add_arrow is TRUE
#'
#' @param y_coord_adjust Vertical adjustment of the arrow origin, default is the same as x_coord_adjust
#'                       Used only if add_arrow is TRUE
#'
#' @param arrow_size The size of the arrow head, default is 0.25 cm
#'                   Used only if add_arrow is TRUE
#'
#' @param arrow_width The width of the arrow, default is 0.2 cm
#'                   Used only if add_arrow is TRUE
#'
#' @param arrow_line_width The width of the arrow line, default is 3
#'                   Used only if add_arrow is TRUE
#'
#' @param arrow_color The color of the arrow and label, default is "black"
#'                    Used only if add_arrow is TRUE
#'
#' @param text_font_size The font size of the text, default is 12
#'                       Used only if add_arrow is TRUE
#'
#' @param text_dst The distance of the text from the arrow in cm, default is 0.3
#'                 Used only if add_arrow is TRUE
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
                              byrow = TRUE,
                              add_arrow = TRUE,
                              bottom_left_only = FALSE,
                              remove_axes = TRUE,
                              remove_legend = FALSE,
                              remove_grid = FALSE,
                              restore_legend = FALSE,
                              legend_position = NULL,
                              legend_title = NULL,
                              legend_title_size = 12,
                              legend_title_bold = FALSE,
                              reduction = "umap",
                              dims = NULL,
                              arrow_length = 2,
                              x_coord_adjust = 0.25,
                              y_coord_adjust = x_coord_adjust,
                              arrow_size = 0.25,
                              arrow_width = 0.2,
                              arrow_line_width = 3,
                              arrow_color = "black",
                              text_font_size = 12,
                              text_dst = 0.3) {
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
    plot <- manipulate_sc_plot(
      plot = plot,
      remove_axes = remove_axes,
      remove_legend = remove_legend,
      remove_grid = remove_grid,
      restore_legend = restore_legend,
      legend_position = legend_position,
      legend_title = legend_title,
      legend_title_size = legend_title_size,
      legend_title_bold = legend_title_bold
    )
    return(plot)
  })

  # Add arrows to the plots
  if (add_arrow && !bottom_left_only) {
    # To be done in all plots
    plot_list <- lapply(plot_list, function(plot) {
      plot <- add_dimension_arrows(
        plot = plot,
        reduction = reduction,
        dims = dims,
        arrow_length = arrow_length,
        x_coord_adjust = x_coord_adjust,
        y_coord_adjust = y_coord_adjust,
        arrow_size = arrow_size,
        arrow_width = arrow_width,
        arrow_line_width = arrow_line_width,
        arrow_color = arrow_color,
        text_font_size = text_font_size,
        text_dst = text_dst
      )
      return(plot)
    })
  } else if (add_arrow && bottom_left_only) {
    # Need to identify which plot is at the bottom left corner
    # Requires: ncol, nrow, byrow
    n_plots <- length(plot_list)
    # If ncol and nrow are not set, then check if the plot is a patchwork object
    if (is.null(nrow) && is.null(ncol)) {
      if (inherits(plots, "patchwork")) {
        nrow <- plots$patches$layout$nrow
        ncol <- plots$patches$layout$ncol
      } else {
        warning("nrow and ncol are not set, and no patchwork default layout is found. Using default ncol = 1 and nrow = length(plots)")
      }
    }

    if (is.null(byrow)) {
      if (inherits(plots, "patchwork") && !is.null(plots$patches$layout$byrow)) {
        byrow <- plots$patches$layout$byrow
      } else {
        byrow <- TRUE
      }
    } else {
      # Do nothing
    }

    # Check which plot is at the bottom left corner
    if (!is.null(ncol)) {
      # Prioritize ncol if set
      nrow <- ceiling(n_plots / ncol)
      if (inherits(plots, "patchwork")) {
        plots$patches$layout$ncol <- ncol
        plots$patches$layout$nrow <- nrow
      }
    } else if (!is.null(nrow)) {
      # Use nrow when ncol is not set
      ncol <- ceiling(n_plots / nrow)
      if (inherits(plots, "patchwork")) {
        plots$patches$layout$nrow <- nrow
        plots$patches$layout$ncol <- ncol
      }
    } else {
      stop("Either ncol or nrow must be set")
    }

    if (byrow) {
      # Fill the plot grid by row
      bottom_left_idx <- (nrow - 1) * ncol + 1
    } else {
      # Fill the plot grid by column
      bottom_left_idx <- nrow
    }
    if (inherits(plots, "patchwork")) {
      plots$patches$layout$byrow <- byrow
    }

    # Add the arrow to the bottom left plot
    plot_list[[bottom_left_idx]] <- add_dimension_arrows(
      plot = plot_list[[bottom_left_idx]],
      reduction = reduction,
      dims = dims,
      arrow_length = arrow_length,
      x_coord_adjust = x_coord_adjust,
      y_coord_adjust = y_coord_adjust,
      arrow_size = arrow_size,
      arrow_width = arrow_width,
      arrow_line_width = arrow_line_width,
      arrow_color = arrow_color,
      text_font_size = text_font_size,
      text_dst = text_dst
    )

  } else {
    # Do nothing
  }

  # Combine the plots into a patchwork object
  if (inherits(plots, "patchwork")) {
    combined_plot <- wrap_plots(plotlist = plot_list,
                                ncol = if (is.null(ncol)) NULL else ncol,
                                nrow = if (is.null(nrow)) NULL else nrow,
                                byrow = byrow,
                                width = plots$patches$layout$width,
                                height = plots$patches$layout$height,
                                guides = plots$patches$layout$guides,
                                tag_level = plots$patches$layout$tag_level,
                                design = plots$patches$layout$design,
                                axes = plots$patches$layout$axes,
                                axis_titles = plots$patches$layout$axis_titles)
    combined_plot <- combined_plot + plots$patches$annotation

  } else if (!is.null(ncol)) {
    combined_plot <- patchwork::wrap_plots(plotlist = plot_list, ncol = ncol, byrow = byrow)
  } else if (!is.null(nrow)) {
    combined_plot <- patchwork::wrap_plots(plotlist = plot_list, nrow = nrow, byrow = byrow)
  } else {
    combined_plot <- patchwork::wrap_plots(plotlist = plot_list)
  }

  return(combined_plot)
}



#'@title Bar Plot for Grouped Data
#'
#' @param df A data frame with columns: `split.by`, `group.by`, and `frequency`.
#'
#' @param group.by Column name used for the fill (segments within bars).
#'
#' @param split.by Column name used for the x-axis (categories to split the bars).
#'
#' @param frequency Column name with the values for bar heights (usually counts
#'                  or proportions).
#'
#' @param mode Either "stack" or "fill".
#'
#' @param palette A qualitative palette name from \code{\link{list_palettes}}
#'                (default "polychrome-36"), or a vector of colors with at least
#'                as many colors as groups.
#'
#' @param title Plot title.
#'
#' @param show_totals Logical, whether to show totals on top of the bars.
#'
#' @return A ggplot2 object
#'
#' @import ggplot2 dplyr
#'
#' @export
#'
grouped_bar_plot <- function(df,
                             group.by,
                             split.by,
                             frequency,
                             mode = c("stack", "fill"),
                             palette = "polychrome-36",
                             title = NULL,
                             show_totals = FALSE) {
  # Check if the required columns are present
  if (!all(c(group.by, split.by, frequency) %in% colnames(df))) {
    stop("Data frame must contain the specified columns.")
  }

  # Check if mode is valid
  mode <- match.arg(mode)

  # Define the colors, from the palette registry or as given
  groups <- df[[group.by]]
  colors <- group_colours(
    palette, if (is.factor(groups)) droplevels(groups) else factor(groups), arg = "palette"
  )$cols

  # Create the bar plot
  plot <- ggplot(df, aes(x = .data[[split.by]], y = .data[[frequency]], fill = .data[[group.by]])) +
    geom_bar(stat = "identity", position = mode) +
    scale_fill_manual(values = colors) +
    theme_classic() +
    labs(y = ifelse(mode == "fill", "Proportion", "Frequency"),
         fill = group.by,
         title = title) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_text(face = "bold"),
      axis.title.x = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  if (show_totals) {
    # Calculate totals for each split.by category
    total_df <- aggregate(df[[frequency]], by = list(df[[split.by]]), FUN = sum)
    colnames(total_df) <- c(split.by, "total")

    # Add totals to the plot
    plot <- plot +
      geom_text(data = total_df,
        aes(x = .data[[split.by]],
                   y = .data[["total"]],
                    label = .data[["total"]]),
        vjust = -0.5,
        fontface = "bold",
        inherit.aes = FALSE,
        size = 3.5
      )
  }

  # Return the plot
  return(plot)
}
