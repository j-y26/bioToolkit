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


#' @title Calculate group proportions and counts across conditions
#'
#' @description Calculates the proportion or count of different groups across
#' specified conditions from a Seurat object or metadata dataframe
#'
#' @param object A Seurat object or a data.frame containing the metadata
#'               (e.g., cell type annotations). If a dataframe is provided,
#'               it must be in long format with columns for each variable, e.g.,
#'               seurat_obj@meta.data.
#'
#' @param group.by Character vector specifying the column name(s) containing the
#'                 grouping variable(s) (e.g., cell types, clusters).
#'
#' @param split.by Character string specifying the column name containing the
#'                 condition variable (e.g., treatment, genotype).
#'
#' @param cells.use Optional character vector specifying the cells to include in the
#'                  analysis. If NULL (default), all cells are included. Can be
#'                  `WhichCells(seurat_obj, idents = c("ident1", "ident2"))` or
#'                  any other vector of cell names.
#'
#' @return A data.frame with columns:
#'   \item{group}{The grouping variable values (factor)}
#'   \item{condition}{The condition/split variable values (factor)}
#'   \item{prop}{The calculated proportion}
#'   \item{count}{The raw count for each group-condition combination (numeric)}
#'   \item{total_count_condition}{The total count for that condition (numeric)}
#'
#' @import dplyr
#'
#' @export
#'
summarize_group_meta <- function(object,
                                 group.by,
                                 split.by,
                                 cells.use = NULL) {
  # Check if the object is a Seurat object or a data frame
  if (is(object, "Seurat")) {
    object <- object@meta.data
  } else if (is.data.frame(object) || inherits(object, "tbl_df")) {
    object <- as.data.frame(object)
  } else {
    stop("The object must be a Seurat object or a data frame-like object")
  }

  # Check if the group.by and split.by columns are in the object
  if (!all(c(group.by, split.by) %in% colnames(object))) {
    stop("The group.by and split.by must be columns in the metadata")
  }

  # Select only the relevant columns
  df <- object[, c(group.by, split.by)]

  # Subset to the specified cells if provided
  if (!is.null(cells.use)) {
    org.ncell <- nrow(df)
    if (!all(cells.use %in% rownames(df))) {
      no_id_cells <- setdiff(cells.use, rownames(df))
      warning(paste0(length(no_id_cells), " cells in cells.use are not found in the object and will be ignored: ", paste(head(no_id_cells, 10), collapse = ", "), if (length(no_id_cells) > 10) "..." else ""))
      df <- df[rownames(df) %in% cells.use, ]
      message(paste0("Subsetting to ", nrow(df), " cells (", round(100 * nrow(df) / org.ncell, 2), "% of original)"))
    } else {
      df <- df[rownames(df) %in% cells.use, ]
      message(paste0("Subsetting to ", nrow(df), " cells (", round(100 * nrow(df) / org.ncell, 2), "% of original)"))
    }
  }

  rm(object) # Free memory
  gc()

  # Rename columns for easier handling
  colnames(df) <- c("group", "condition")

  # Process NA values
  if (any(is.na(df))) {
    na_idx <- is.na(df$group) | is.na(df$condition)
    warning(paste0(sum(na_idx), " cells with NA values in group.by or split.by will be ignored"))
    df <- df[!na_idx, ]
    message(paste0("Remaining ", nrow(df), " cells after removing NAs"))
  }

  # Ensure there are still cells left
  if (nrow(df) == 0) {
    stop("No cells remaining after filtering. Please check your inputs.")
  }

  # Calculate counts
  summary_df <- table(df$group, df$condition) %>%
    as.data.frame() %>%
    rename(group = Var1, condition = Var2, count = Freq)

  # Calculate total counts per condition
  total_counts <- summary_df %>%
    group_by(condition) %>%
    summarise(total_count_condition = sum(count))

  # Calculate proportions
  summary_df <- summary_df %>%
    left_join(total_counts, by = "condition") %>%
    mutate(prop = count / total_count_condition)

  # Re-arrange columns
  summary_df <- summary_df %>%
    select(group, condition, prop, count, total_count_condition)

  return(summary_df)
}


#' @title Metadata Scatter Plot
#'
#' @description This function generates a scatter plot for visualizing the
#'              relationship of a grouped metadata variable between two conditions.
#'              Each point represents a group (e.g., cell type) and its position
#'              is determined by the proportion or count of that group in each
#'              of two specified conditions. This is useful for comparing how
#'              the distribution of groups changes between conditions.
#'
#' @param object A Seurat object or a data.frame containing the metadata
#'               (e.g., cell type annotations). If a dataframe is provided,
#'               it must be in long format with columns for each variable, e.g.,
#'               `seurat_obj@meta.data`.
#'
#' @param group.by Character vector specifying the column name(s) containing the
#'                 grouping variable(s) (e.g., cell types, clusters).
#'
#' @param split.by Character string specifying the column name containing the
#'                 condition variable (e.g., treatment, genotype).
#'
#' @param cells.use Optional character vector specifying the cells to include in the
#'                  analysis. If NULL (default), all cells are included. Can be
#'                  `WhichCells(seurat_obj, idents = c("ident1", "ident2"))` or
#'                  any other vector of cell names.
#'
#' @param group.use Optional character vector specifying the groups to include
#'                  in the plot. If NULL (default), all groups are included.
#'                  When plotting proportions, setting this parameter does NOT
#'                  affect the proportion calculations, which are always based
#'                  on all groups present in the data. This will only filter the
#'                  points shown in the plot.
#'
#' @param x_condition The name of the condition to plot on the x-axis. Must be
#'                    one of the unique values in the split.by column. Optional
#'                   (default NULL) if there are only two unique conditions.
#'
#' @param y_condition The name of the condition to plot on the y-axis. Must be
#'                    one of the unique values in the split.by column. Optional
#'                   (default NULL) if there are only two unique conditions.
#'
#' @param x_label The label for the x-axis. If NULL (default), the x_condition
#'                 will be used as the label.
#'
#' @param y_label The label for the y-axis. If NULL (default), the y_condition
#'                 will be used as the label.
#'
#' @param mode Either "prop" or "count" to determine whether to plot the
#'             proportion or count of cells for each group in the specified
#'             conditions.
#'
#' @param percentage Logical, whether to convert proportions to percentages
#'                   (default FALSE).
#'
#' @param colors A named vector of colors to use for the groups. The names
#'              should match the group names and should have exactly the same
#'              length as the number of unique groups.
#'
#' @param point_size The size of the points in the scatter plot, default is 3.
#'
#' @param point_alpha The transparency of the points in the scatter plot,
#'                   default is 0.8.
#'
#' @param point_border Whether to add a black border around the points,
#'                     default is TRUE.
#'
#' @param point_border_stroke The thickness of the point border,
#'                             default is 0.5.
#'
#' @param add_diagonal Logical, whether to add a diagonal reference line
#'                  (default TRUE).
#'
#' @param label Logical, whether to add labels to the points (default TRUE).
#'
#' @param label_text_size The size of the point labels, default is 3.5.
#'
#' @param legend_position The position of the legend, default is 'none'.
#'
#' @param legend_title The title of the legend, default is the value of group.by.
#'
#' @param remove_zero Logical, whether to remove groups with zero counts in both
#'                    conditions (default TRUE).
#'
#'
#' @return A ggplot object representing the scatter plot.
#'
#' @import ggplot2 dplyr
#'
#' @export
#'
plot_meta_scatter <- function(
  object,
  group.by,
  split.by,
  cells.use = NULL,
  group.use = NULL,
  x_condition = NULL,
  y_condition = NULL,
  x_label = NULL,
  y_label = NULL,
  mode = c("prop", "count"),
  percentage = TRUE,
  colors = NULL,
  point_size = 4,
  point_alpha = 0.8,
  point_border = TRUE,
  point_border_stroke = 1,
  add_diagonal = TRUE,
  diagonal_slope = 1,
  diagonal_intercept = 0,
  diagonal_color = "grey60",
  diagonal_linetype = "11",
  diagonal_size = 0.8,
  label = TRUE,
  label_text_size = 4,
  legend_position = "none",
  legend_title = group.by,
  remove_zero = TRUE
) {
  # Get the calculated summary data
  summary_df <- summarize_group_meta(object, group.by, split.by, cells.use)

  # Pull metadata if a Seurat object was provided
  if (is(object, "Seurat")) {
    object <- object@meta.data
  } else if (is.data.frame(object) || inherits(object, "tbl_df")) {
    object <- as.data.frame(object)
  } else {
    stop("The object must be a Seurat object or a data frame-like object")
  }

  # Mode
  mode <- mode[1]
  if (!mode %in% c("prop", "count")) {
    stop("Mode must be either 'prop' or 'count'.")
  }

  # Check if the specified conditions are valid
  if (is.null(levels(object[[split.by]]))) {
    conditions <- unique(object[[split.by]])
  } else {
    conditions <- levels(object[[split.by]])
  }
  # If conditions is a dataframe, convert to character vector
  if (is.data.frame(conditions)) {
    conditions <- conditions[[split.by]]
  }
  unique_conditions <- conditions[conditions %in% unique(summary_df$condition)] %>% unique()
  if (length(unique_conditions) < 2) {
    stop("The split.by column must contain at least two unique conditions.")
  }
  if (is.null(x_condition) && is.null(y_condition)) {
    if (length(unique_conditions) > 2) {
      stop("When there are more than two unique conditions, both x_condition and y_condition must be specified.")
    } else {
      x_condition <- unique_conditions[1]
      y_condition <- unique_conditions[2]
      message(paste0("Using ", x_condition, " for x-axis and ", y_condition, " for y-axis."))
    }
  } else if (is.null(x_condition) || is.null(y_condition)) {
    stop("Both x_condition and y_condition must be specified or both left as NULL.")
  } else {
    if (!(x_condition %in% unique_conditions)) {
      stop("x_condition must be one of the unique values in the split.by column.")
    }
    if (!(y_condition %in% unique_conditions)) {
      stop("y_condition must be one of the unique values in the split.by column.")
    }
    if (x_condition == y_condition) {
      stop("x_condition and y_condition must be different.")
    }
  }

  # Retain only the specified conditions
  summary_df <- summary_df %>% filter(condition %in% c(x_condition, y_condition))
  summary_df$condition <- factor(summary_df$condition, levels = c(x_condition, y_condition))

  # Check if the specified groups are valid
  if (!is.null(group.use)) {
    if (!all(group.use %in% unique(summary_df$group))) {
      stop("All values in group.use must be present in the group.by column.")
    }
    summary_df <- summary_df %>% filter(group %in% group.use)
  }

  if (is.null(levels(object[[group.by]]))) {
    groups <- unique(object[[group.by]])
  } else {
    groups <- levels(object[[group.by]])
  }
  # If groups is a dataframe, convert to character vector
  if (is.data.frame(groups)) {
    groups <- groups[[group.by]]
  }
  unique_groups <- groups[groups %in% unique(summary_df$group)] %>% unique()
  if (length(unique_groups) < 1) {
    stop("The group.by column must contain at least one valid group.")
  }
  summary_df$group <- factor(summary_df$group, levels = unique_groups)

  # Get the x and y data for plotting
  if (mode == "prop") {
    x_data <- summary_df %>%
        filter(condition == x_condition) %>%
        filter(group %in% unique_groups) %>%
        select(group, prop) %>%
        rename(x_value = prop)
    y_data <- summary_df %>%
        filter(condition == y_condition) %>%
        filter(group %in% unique_groups) %>%
        select(group, prop) %>%
        rename(y_value = prop)
    if  (percentage) {
      x_data$x_value <- x_data$x_value * 100
      y_data$y_value <- y_data$y_value * 100
    }
  } else {
    x_data <- summary_df %>% filter(condition == x_condition) %>% select(group, count) %>% rename(x_value = count)
    y_data <- summary_df %>% filter(condition == y_condition) %>% select(group, count) %>% rename(y_value = count)
  }

  plot_df <- merge(x_data, y_data, by = "group")
  plot_df[is.na(plot_df)] <- 0 # Replace NA with 0
  rownames(plot_df) <- plot_df$group
  plot_df <- plot_df[, c("x_value", "y_value", "group")]

  # Define axis labels
  if (is.null(x_label)) {
    x_label <- paste0(if (mode == "prop" && percentage) "Percentage" else if (mode == "prop" && !percentage) "Proportion" else "Count", " in ", x_condition)
  }
  if (is.null(y_label)) {
    y_label <- paste0(if (mode == "prop" && percentage) "Percentage" else if (mode == "prop" && !percentage) "Proportion" else "Count", " in ", y_condition)
  }

  # Remove groups with zero counts in both conditions if requested
  if (remove_zero) {
    zero_idx <- plot_df$x_value == 0 & plot_df$y_value == 0
    if (any(zero_idx)) {
      warning(paste0(sum(zero_idx), " group(s) with zero counts in both conditions will be removed."))
      plot_df <- plot_df[!zero_idx, ]
      message(paste0("Remaining ", nrow(plot_df), " group(s) after removing zeros."))

      # Update factor levels
      plot_df$group <- factor(plot_df$group)
      # Update colors if provided
      if (!is.null(colors)) {
        colors <- colors[levels(plot_df$group)]
      }

      if (nrow(plot_df) == 0) {
        stop("No groups remaining after removing zeros. Please check your inputs.")
      }
    }
  }

  rm(object) # Free memory
  gc()

  # Generate the base plot
  p <- ggplot(plot_df, aes(x = x_value, y = y_value, fill = group)) +
    geom_point(
      size = point_size,
      alpha = point_alpha,
      shape = ifelse(point_border, 21, 16),
      stroke = ifelse(point_border, point_border_stroke, 0),
      color = ifelse(point_border, "black", NA)
    ) +
    theme_classic() +
    labs(
      x = x_label,
      y = y_label
    ) +
    theme(
      axis.title = element_text(face = "bold"),
      legend.position = legend_position
    )

  # Check legend
  if (legend_position != "none") {
    p <- p + theme(
      legend.title = element_text(face = "bold")
    ) +
      labs(fill = legend_title)
  }

  # Add diagonal reference line if requested
  if (add_diagonal) {
    p <- p + geom_abline(slope = diagonal_slope, intercept = diagonal_intercept, linetype = diagonal_linetype, color = diagonal_color, size = diagonal_size)
  }

  # Add custom colors if provided
  if (!is.null(colors)) {
    if (length(colors) != length(unique(plot_df$group))) {
      stop("The colors must have the same length as the number of unique groups.")
    } else if (!all(names(colors) %in% unique(plot_df$group))) {
      stop("The colors must have the same names as the unique groups.")
    } else {
      # arrange the colors based on the groups
      colors <- colors[unique(plot_df$group)]
      p <- p + scale_fill_manual(values = colors)
    }
  }

  # Add labels if requested
  if (label) {
    check_installed("ggrepel", reason = "for repelled text labels")
    p <- p + ggrepel::geom_text_repel(
      aes(label = group),
      color = "black",
      size = label_text_size,
      segment.color = "grey50",
      segment.size = 0.6,
      max.overlaps = Inf)
  }

  return(p)
}
