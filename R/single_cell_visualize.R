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
#' @param down_color The color for down-regulated features, default is "darkgreen"
#'
#' @param mid_color The color for 0 logFC features, default is "white"
#'
#' @param up_color The color for up-regulated features, default is "darkslateblue"
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
  down_color = "darkgreen",
  mid_color = "white",
  up_color = "darkslateblue",
  nsig_color = "grey70",
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
  x_group_label_position = "bottom",
  y_group_label_position = "left"

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
        color = ifelse(signif, !!sym(logFC_var), NA)
      ))
  } else {
    dot_plot <- ggplot(de_results,
      aes(x = !!sym(comparison_var),
          y = !!sym(feature_var))) +
      geom_point(aes(
        size = dot_size,
        color = ifelse(signif, !!sym(logFC_var), NA)
      ))
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
    scale_color_gradient2(
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
