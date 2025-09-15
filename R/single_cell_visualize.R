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

  # Get the dimension names if not provided
  if (is.null(dims)) {
    # Check if the reduction method is valid
      reduction <- tolower(reduction)
    if (grepl("umap", reduction)) {
      dims <- c("UMAP-1", "UMAP-2")
    } else if (grepl("tsne", reduction) | grepl("t-sne", reduction)){
      dims <- c("tSNE-1", "tSNE-2")
    } else if (grepl("pca", reduction)) {
      dims <- c("PC-1", "PC-2")
    } else if (grepl("spatial", reduction)) {
      dims <- c("Spatial-1", "Spatial-2")
    } else {
      dims <- c(paste0(reduction, "-1"), paste0(reduction, "-2"))
    }
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
  arrow_grob <- grobTree(arrow_x, arrow_y, arrowhead_x, arrowhead_y, label_x, label_y)

  plot <- plot +
    annotation_custom(arrow_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    coord_cartesian(clip = "off") +
    theme(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank()
    )

  return(plot)
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
#' @import Seurat
#'
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
    plot <- plot + NoAxes()
  }
  # Remove legend
  if (remove_legend) {
    plot <- plot + NoLegend()
  }
  # Remove grid
  if (remove_grid) {
    plot <- plot + NoGrid()
  }
  # Restore legend
  if (restore_legend) {
    plot <- plot + RestoreLegend()
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
#' @import Seurat
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


#' @title Dot Plot for comparisons with logFC and p-values
#'
#' @description This function generates a dot plot for visualizing comparisons
#'              between multiple pairs of groups. The dot plot shows the log
#'              fold change (logFC) and p-values for each comparison for each
#'              feature. The dot plot is useful for visualizing differential
#'              expression analysis results. If a gene in a given comparison has
#'              a significant p-value, the dot is colored based on the logFC
#'              value. The size of the dot is proportional to the -log10 of the
#'              p-value. If not significant, the dot is colored in grey.
#'
#' @param de_results A data frame with the differential expression results, must
#'                   contain the columns for "feature", "avg_logFC", "p_val_adj",
#'                   and "comparison"
#'
#' @param features A vector of feature names to plot
#'
#' @param feature_var The name of the feature variable in the de_results data
#'                    frame, default is "gene"
#'
#' @param logFC_var The name of the logFC variable in the de_results data frame,
#'                  default is "avg_logFC"
#'
#' @param p_val_var The name of the p-value variable in the de_results data frame,
#'                  default is "p_val_adj"
#'
#' @param comparison_var The name of the comparison variable in the de_results
#'                       data frame, default is "comparison"
#'
#' @param logFC_cutoff The logFC cutoff for filtering the features, default is 0
#'
#' @param p_val_cutoff The p-value cutoff for filtering the features, default is
#'                     0.05
#'
#' @param palette The color palette to use for the dot plot for representing the
#'                logFC values, default is "NULL". If NULL, down_color, mid_color,
#'                and up_color will be used. If not NULL, the palette will be used
#'                to represent the logFC values and the down_color, mid_color, and
#'                up_color will be ignored
#'
#' @param down_color The color for down-regulated features
#'
#' @param mid_color The color for 0 logFC features
#'
#' @param up_color The color for up-regulated features
#'
#' @param nsig_color The color for non-significant features, default is "grey70"
#'
#' @param col_limits_logFC The limits for the logFC color scale, default is
#'                         using the min and max logFC values in the data.
#'                         This argument is useful for when extreme logFC values
#'                         are present in the data, which could result in
#'                         most dots share very vague colors. The limits can
#'                         be set to a smaller range to enhance the color
#'                         contrast. The limits should be a vector of two
#'                         numeric values, e.g., c(-2, 2)
#'
#' @param sig_only Whether to plot only significant features, default is FALSE,
#'                 which plots all features provided. If TRUE, the plot will
#'                 only show features-comparisons with p-values < p_val_cutoff
#'                 and logFC > logFC_cutoff
#'
#' @param invert Whether to invert the plot, default is FALSE. By default, the
#'               x axis is the features and the y axis is the comparisons. If
#'               TRUE, the x axis is the comparisons and the y axis is the
#'               features
#'
#' @param x_text_angle The angle of the x-axis text, default is 45
#'
#' @param group_by_features The column name in the de_results data frame to
#'                          to specify the grouping of the features, default is
#'                          NULL. Having only unique values in the
#'                          group_by_features is NOT allowed. To handle the
#'                          order in which the groups are plotted, set the
#'                          group_by_features to a factor with the desired order
#'                          using \code{factor(group_by_features, levels = c("levels"))}
#'
#' @param group_by_comparisons The column name in the de_results data frame to
#'                             to specify the grouping of the comparisons,
#'                             default is NULL. Having only unique values in the
#'                             group_by_comparisons is NOT allowed. To handle the
#'                             order in which the groups are plotted, set the
#'                             group_by_comparisons to a factor with the desired
#'                             order using \code{factor(group_by_comparisons, levels = c("levels"))}
#'
#' @param panel_spacing The spacing between the panels, default is 0.2
#'
#' @param legend_position The position of the legend, default is 'right'
#'
#' @param legend_key_size The size of the legend key, default is 0.8
#'
#' @param legend_title The title of the legend, default is "Avg log2(FC)"
#'
#' @param legend_title_size The size of the legend title, default is 8
#'
#' @param legend_text_size The size of the legend text, default is 7
#'
#' @param x_group_label_position The position of the x-axis group label, default is "bottom"
#'
#' @param y_group_label_position The position of the y-axis group label, default is "left"
#'
#' @param dot_border The border size of the dots, default is 0.3
#'
#' @return a ggplot object with the dot plot
#'
#' @import ggplot2 dplyr
#'
#' @importFrom RColorBrewer brewer.pal
#'
#' @importFrom scales gradient_n_pal
#'
#' @export
#'
de_dot_plot <- function(
  de_results,
  features,
  feature_var = "gene",
  logFC_var = "avg_log2FC",
  p_val_var = "p_val_adj",
  comparison_var = "comparison",
  logFC_cutoff = 0,
  p_val_cutoff = 0.05,
  palette = NULL,
  down_color = "#00441B",
  mid_color = "#F0F0F0",
  up_color = "#2A004F",
  nsig_color = "grey30",
  col_limits_logFC = c(min(de_results[[logFC_var]]), max(de_results[[logFC_var]])),
  sig_only = FALSE,
  invert = FALSE,
  x_text_angle = 45,
  group_by_features = NULL,
  group_by_comparisons = NULL,
  panel_spacing = 0.2,
  legend_position = "right",
  legend_key_size = 0.8,
  legend_title = "Avg log2(FC)",
  legend_title_size = 8,
  legend_text_size = 7,
  x_group_label_position = "top",
  y_group_label_position = "left",
  dot_border = 0.3

) {
  # Check if the de_results is a data frame
  if (!is.data.frame(de_results)) {
    stop("The de_results must be a data frame")
  }

  # Check if the feature_var is in the de_results
  if (!feature_var %in% colnames(de_results)) {
    stop(paste0("The feature_var ", feature_var, " must be a column in the de_results data frame"))
  }

  # Check if the logFC_var is in the de_results
  if (!logFC_var %in% colnames(de_results)) {
    stop(paste0("The logFC_var ", logFC_var, " must be a column in the de_results data frame"))
  }

  # Check if the p_val_var is in the de_results
  if (!p_val_var %in% colnames(de_results)) {
    stop(paste0("The p_val_var ", p_val_var, " must be a column in the de_results data frame"))
  }

  # Check if the comparison_var is in the de_results
  if (!comparison_var %in% colnames(de_results)) {
    stop(paste0("The comparison_var ", comparison_var, " must be a column in the de_results data frame"))
  }

  # Check if the logFC_cutoff is numeric
  if (!is.numeric(logFC_cutoff)) {
    stop("The logFC_cutoff must be a numeric value")
  }

  # Check if the p_val_cutoff is numeric
  if (!is.numeric(p_val_cutoff)) {
    stop("The p_val_cutoff must be a numeric value")
  }

  #=============================================================================
  # Currently does not support having BOTH group_by_features and group_by_comparisons
  #=============================================================================
  if (!is.null(group_by_features) && !is.null(group_by_comparisons)) {
    stop("Currently does not support having BOTH group_by_features and group_by_comparisons")
  }

  # Validate grouping columns
  if (!is.null(group_by_features) && !group_by_features %in% colnames(de_results)) {
    stop(paste0("The group_by_features ", group_by_features, " must be a column in the de_results data frame"))

    if (length(unique(de_results[[group_by_features]]) == 1)) {
      stop(paste0("The group_by_features ", group_by_features, " must have more than one unique value"))
    }
  }

  if (!is.null(group_by_comparisons) && !group_by_comparisons %in% colnames(de_results)) {
    stop(paste0("The group_by_comparisons ", group_by_comparisons, " must be a column in the de_results data frame"))

    if (length(unique(de_results[[group_by_comparisons]]) == 1)) {
      stop(paste0("The group_by_comparisons ", group_by_comparisons, " must have more than one unique value"))
    }
  }

  # Find and note any requested features that are not in the de_results
  # Filter the de_results to only include the requested features
  missing_features <- setdiff(features, unique(de_results[[feature_var]]))
  if (length(missing_features) > 0) {
    warning(paste0("The following features are not in the de_results: ", paste(missing_features, collapse = ", ")))
  }
  de_results <- de_results %>% filter(tolower(!!sym(feature_var)) %in% tolower(features))

  if (nrow(de_results) == 0) {
    stop("No matching features found in the provided data frame.")
  }

  # Add a signif column
  de_results$signif <- de_results[[p_val_var]] < p_val_cutoff & abs(de_results[[logFC_var]]) > logFC_cutoff

  # Filter the features based on the logFC and p-value cutoffs
  if (sig_only) {
    de_results <- de_results %>% filter(signif)
  }

  # Alter if any features are removed
  removed_features <- setdiff(features, unique(de_results[[feature_var]]))
  if (length(removed_features) > 0) {
    warning(paste0("The following features were removed due to  non-significance: ", paste(removed_features, collapse = ", ")))
  }

  # Compute the dot size
  de_results$dot_size <- -log10(de_results[[p_val_var]])

# Create the base plot with invert option and ordered axes
  if (!invert) {
    dot_plot <- ggplot(de_results,
      aes(x = !!sym(feature_var),
          y = !!sym(comparison_var))) +
      geom_point(aes(
        size = dot_size,
        fill = !!sym(logFC_var)
      ),
        shape = 21,
        stroke = dot_border,
        color = "black"
      )
  } else {
    dot_plot <- ggplot(de_results,
      aes(x = !!sym(comparison_var),
          y = !!sym(feature_var))) +
      geom_point(aes(
        size = dot_size,
        fill = !!sym(logFC_var)
      ),
        shape = 21,
        stroke = dot_border,
        color = "black"
      )
  }

  # Add faceting based on grouping variables
  if (!is.null(group_by_features)) {
    if (!invert) {
      dot_plot <- dot_plot +
        facet_grid(as.formula(paste(".", "~", group_by_features)),
                    scales = "free_x",
                    space = "free_x",
                    switch = if (x_group_label_position == "top") NULL else "x")
    } else {
      dot_plot <- dot_plot +
        facet_grid(as.formula(paste(group_by_features, "~ .")),
                    scales = "free_y",
                    space = "free_y",
                    switch = if (y_group_label_position == "left") "y" else NULL)
    }
  } else if (!is.null(group_by_comparisons)) {
    if (!invert) {
      dot_plot <- dot_plot +
        facet_grid(as.formula(paste(group_by_comparisons, "~ .")),
                    scales = "free_y",
                    space = "free_y",
                    switch = if (y_group_label_position == "left") "y" else NULL)
    } else {
      dot_plot <- dot_plot +
        facet_grid(as.formula(paste(".", "~", group_by_comparisons)),
                    scales = "free_x",
                    space = "free_x",
                    switch = if (y_group_label_position == "left") "x" else NULL)
    }
  }

  # Check the colors for the dots
  if (!is.null(palette)) {
    # Obtain the low, mid, and high colors from the palette
    palette_colors <- RColorBrewer::brewer.pal(9, palette)
    pal_func <- scales::gradient_n_pal(palette_colors)
    down_color <- pal_func(0)
    mid_color <- pal_func(0.5)
    up_color <- pal_func(1)
  }

  # Add the color scale and point size
  dot_plot <- dot_plot  +
    scale_fill_gradient2(
      low = down_color,
      mid = mid_color,
      high = up_color,
      midpoint = 0,
      na.value = nsig_color,
      limits = col_limits_logFC,
      name = legend_title
    ) +
    scale_size_continuous(name = "-log10(p-value)") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = x_text_angle, hjust = 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", vjust = 1),
      strip.placement = "outside",
      panel.spacing = unit(panel_spacing, "lines"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.key.size = unit(legend_key_size, "lines"),
      legend.margin = margin(0, 0, 0, 0),
      legend.position = legend_position,
      legend.title = element_text(size = legend_title_size),
      legend.text = element_text(size = legend_text_size),
      legend.box.margin = margin(-10, 0, 0, 0),
    )

  return(dot_plot)
}


#' @title Annotation Sankey Plot
#'
#' @description This function generates a Sankey plot for visualizing the
#'              transitions and overlap between different annotations. The
#'              Sankey plot is useful for visualizing the flow of cells
#'              between different levels of annotations.
#'
#' @param object A Seurat object or a data frame (\code{suerat_obj@meta.data})
#'               with the annotations to visualize
#'
#' @param annotations A vector of annotation names to visualize. The order will
#'                    dictate the flow in the Sankey plot
#'
#' @param title The title of the Sankey plot, default is 'Annotation Sankey Plot'
#'
#' @param nodes A vector of unique node/group names to use in the Sankey
#'              plot. If NULL, all unique groups in the annotations will be used
#'
#' @param colors A named vector of colors to use for the annotations. The names
#'               should match the annotation group names and should have exactly
#'               the same length as the unique groups in the annotations
#'
#' @param thickness The thickness of the Sankey plot, default is 15
#'
#' @param pad The padding between the nodes, default is 15
#'
#' @param font_size The font size of the text, default is 10
#'
#' @importFrom plotly plot_ly layout
#' @import dplyr
#'
#' @export
#'
#' @return a plotly object with the Sankey plot
#'
annotation_sankey_plot <- function(
  object,
  annotations,
  title = "Annotation Sankey Plot",
  nodes = NULL,
  colors = NULL,
  thickness = 15,
  pad = 15,
  font_size = 10
) {
  # Check if the object is a Seurat object or a data frame
  if (is(object, "Seurat")) {
    object <- object@meta.data
  } else if (is.data.frame(object)) {
    object <- object
  } else {
    stop("The object must be a Seurat object or a data frame")
  }

  # Check annotation columns
  if (length(annotations) < 2) {
    stop("The annotations must contain at least two columns")
  }

  # Check if the annotations are in the object
  if (!all(annotations %in% colnames(object))) {
    stop("The annotations must be columns in the object")
  }

  # Filter to keep only the annotations
  object <- object[, c(annotations)]

  # Filter nodes
  if (is.null(nodes)) {
    nodes <- unique(unlist(object))
  } else {
    object <- object %>% filter_all(any_vars(. %in% nodes))
  }

  # Assign colors to the nodes
  if (is.null(colors)) {
    colors <- scales::hue_pal()(length(nodes))
  } else {
    if (length(colors) != length(nodes)) {
      stop("The colors must have the same length as the nodes")
    } else if (!all(names(colors) %in% nodes)) {
      stop("The colors must have the same names as the nodes")
    } else {
      # arrange the colors based on the nodes
      colors <- colors[nodes]
    }
  }

  # Generate list that represents the links
  links <- list()
  for (i in 1:(length(annotations) - 1)) {
    # Get the unique groups in the current and next annotation
    groups_current <- unique(object[[annotations[i]]])
    groups_next <- unique(object[[annotations[i + 1]]])

    # Generate the links
    for (group_current in groups_current) {
      for (group_next in groups_next) {
        link <- list(
          source = group_current,
          target = group_next,
          value = sum((object[[annotations[i]]] == group_current) & (object[[annotations[i + 1]]] == group_next))
        )
        if (!link$value == 0) {
          links$source <- c(links$source, link$source)
          links$target <- c(links$target, link$target)
          links$value <- c(links$value, link$value)
        }
      }
    }
  }

  # Generate the sankey plot
  plt <- plotly::plot_ly(
    type = "sankey",
    domain = list(
      x =  c(0,1),
      y =  c(0,1)
    ),
    orientation = "h",
    valueformat = ".0f",
    valuesuffix = "TWh",
    node = list(
      label = nodes,
      color = colors,
      pad = pad,
      thickness = thickness,
      line = list(
        color = "black",
        width = 0.5
      )
    ),

    link = links
  )

  # Add layout
  plt <- plt %>% plotly::layout(
    title = title,
    font = list(
      size = font_size
    ),
    xaxis = list(
      showgrid = FALSE,
      zeroline = FALSE
    ),
    yaxis = list(
      showgrid = FALSE,
      zeroline = FALSE
    )
  )

  return(plt)
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
#' @param palette A valid palette name from "alphabet", "alphabet2", "glasbey", "polychrome", "stepped", and "parade".
#'                Alternatively, a vector of colors can be provided with the same length as unique groups.
#'
#' @param title Plot title.
#'
#' @param show_totals Logical, whether to show totals on top of the bars.
#'
#' @return A ggplot2 object
#'
#' @import ggplot2 dplyr
#' @importFrom Seurat DiscretePalette
#'
#' @export
#'
grouped_bar_plot <- function(df,
                             group.by,
                             split.by,
                             frequency,
                             mode = c("stack", "fill"),
                             palette = "polychrome",
                             title = NULL,
                             show_totals = FALSE) {
  # Check if the required columns are present
  if (!all(c(group.by, split.by, frequency) %in% colnames(df))) {
    stop("Data frame must contain the specified columns.")
  }

  # Check if mode is valid
  mode <- match.arg(mode)
  if (!mode %in% c("stack", "fill")) {
    stop("Mode must be either 'stack' or 'fill'.")
  }

  # Define the colors
  if (length(palette) == 1 && palette %in% c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped", "parade")) {
    colors <- Seurat::DiscretePalette(
      n = length(unique(df[[group.by]])),
      palette = palette
    )
  } else if (length(palette) == length(unique(df[[group.by]]))) {
    # Use custom colors
    colors <- palette
  } else {
    stop("Palette must be one of 'alphabet', 'alphabet2', 'glasbey', 'polychrome', 'stepped', and 'parade' or a vector of colors with the same length as unique groups.")
  }

  # Create the bar plot
  plot <- ggplot(df, aes_string(x = split.by, y = frequency, fill = group.by)) +
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
        aes_string(x = split.by,
                   y = "total",
                    label = "total"),
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
summarize_group_meta <- function(object, group.by, split.by, cells.use = NULL) {
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

  rm(object) # Free memory
  gc()

  # Subset to the specified cells if provided
  if (!is.null(cells.use)) {
    org.ncell <- nrow(df)
    if (!all(cells.use %in% rownames(object))) {
      no_id_cells <- setdiff(cells.use, rownames(object))
      warning(paste0(length(no_id_cells), " cells in cells.use are not found in the object and will be ignored: ", paste(head(no_id_cells, 10), collapse = ", "), if (length(no_id_cells) > 10) "..." else ""))
      df <- df[rownames(df) %in% cells.use, ]
      message(paste0("Subsetting to ", nrow(df), " cells (", round(100 * nrow(df) / org.ncell, 2), "% of original)"))
    } else {
      df <- df[rownames(df) %in% cells.use, ]
      message(paste0("Subsetting to ", nrow(df), " cells (", round(100 * nrow(df) / org.ncell, 2), "% of original)"))
    }
  }

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
#' @importFrom ggrepel geom_text_repel
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
