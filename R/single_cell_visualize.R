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
#' @param arrow_length The size of the arrow, default is 0.15 (15% of the plot)
#'
#' @param arrow_x_adjust_ratio The x-axis adjustment ratio for the arrow,
#'                             increase (> 1) and decrease (< 1) the x-axis.
#'                             Default is 1
#'
#' @param arrow_y_adjust_ratio The y-axis adjustment ratio for the arrow,
#'                             increase (> 1) and decrease (< 1) the y-axis.
#'                             Default is 1
#'
#' @param x_coord_adjust Horizontal adjustment of the arrow origin, default is 0
#'
#' @param y_coord_adjust Vertical adjustment of the arrow origin, default is 0
#'
#' @param arrow_size The size of the arrow head, default is 0.2 cm
#'
#' @param arrow_width The width of the arrow, default is 0.9
#'
#' @param arrow_text_color The color of the arrow and label, default is "black"
#'
#' @param arrow_type The type of the arrow, default is "closed"
#'
#' @param text_font_size The font size of the text, default is 4
#'
#' @param relative_text_dst The relative distance of the text from the arrow,
#'                          default is 0.15
#'
#'
#' @return a ggplot object with the dimension arrows added
#'
#'
#' @import ggplot2
#'
#' @export
#'
add_dimension_arrows <- function(plot,
                                 reduction = "umap",
                                 dims = NULL,
                                 arrow_length = 0.1,
                                 arrow_x_adjust_ratio = 1,
                                 arrow_y_adjust_ratio = 1,
                                 x_coord_adjust = 0,
                                 y_coord_adjust = 0,
                                 arrow_size = 0.2,
                                 arrow_width = 0.9,
                                 arrow_text_color = "black",
                                 arrow_type = "closed",
                                 text_font_size = 4,
                                 relative_text_dst = 0.15) {
  # Check if the plot is a ggplot object
  if (!is.ggplot(plot)) {
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
      stop("Invalid reduction method, allowed values must contain 'umap', 'tsne', 'pca', 'spatial'")
    }
  }

  # Extract plot data and compute arrow length
  plot_data <- ggplot_build(plot)
  x_range <- plot_data$layout$panel_scales_x[[1]]$range$range
  y_range <- plot_data$layout$panel_scales_y[[1]]$range$range

  # Compute the base arrow length
  arrow_length <- arrow_size * min(diff(x_range), diff(y_range))

  # Compute the arrow coordinates
  arrow_x <- x_range[1] + x_coord_adjust
  arrow_y <- y_range[1] + y_coord_adjust

  # Add the arrows to the plot
  plot <- plot +
    geom_segment(
      x = arrow_x,
      y = arrow_y,
      xend = arrow_x + arrow_length * arrow_x_adjust_ratio,
      yend = arrow_y,
      arrow = arrow(type = arrow_type, length = unit(arrow_size, "cm")),
      lwd = arrow_width,
      color = arrow_text_color
    ) +
    geom_segment(
      x = arrow_x,
      y = arrow_y,
      xend = arrow_x,
      yend = arrow_y + arrow_length * arrow_y_adjust_ratio,
      arrow = arrow(type = arrow_type, length = unit(arrow_size, "cm")),
      lwd = arrow_width,
      color = arrow_text_color
    ) +
    geom_text(
      x = arrow_x + (arrow_length * arrow_x_adjust_ratio) / 2,
      y = arrow_y - arrow_length * relative_text_dst,
      label = dims[1],
      size = text_font_size,
      color = arrow_text_color
    ) +
    geom_text(
      x = arrow_x - arrow_length * relative_text_dst,
      y = arrow_y + (arrow_length * arrow_y_adjust_ratio) / 2,
      label = dims[2],
      size = text_font_size,
      color = arrow_text_color,
      angle = 90
    )

  return(plot)
}



#' @title do_DimPlot with arrows on the bottom left corner
#'
#' @description This function is a further implementation of the do_DimPlot
#'              function from the SCpubr package. It generates a ggplot object
#'              with the dimension arrows added to the bottom left corner of
#'              the plot. The function is useful for visualizing single cell
#'              data with dimension reduction methods such as UMAP, t-SNE, PCA,
#'              and spatial coordinates.
#'
#' @param object A Seurat object
#'
#' @param reduction The reduction method used to generate the plot, default is
#'                 "umap".
#'
#' @param dims A vector of dimension names, overwrites the default dimension
#'             names assumed from reduction method
#'
#' @param group.by A vector of cell metadata to group the cells by, default is
#'                 NULL
#'
#' @param colors.use Named vector of colors to use for plotting, default is NULL.
#'
#' @param pt.size The size of the points, default is 0.5
#'
#' @param shuffle Whether to shuffle the dots. Default is FALSE to be consistent
#'                with Seurat
#'
#' @param label Whether to label the groups, default is FALSE
#'
#' @param repel Whether to repel the labels, default is to be consistent with
#'              the "label" argument
#'
#' @param legend.position The position of the legend, default is 'right'
#'
#' @param arrow_length The size of the arrow, default is 0.15 (15% of the plot)
#'
#' @param arrow_x_adjust_ratio The x-axis adjustment ratio for the arrow,
#'                             increase (> 1) and decrease (< 1) the x-axis.
#'                             Default is 1
#'
#' @param arrow_y_adjust_ratio The y-axis adjustment ratio for the arrow,
#'                             increase (> 1) and decrease (< 1) the y-axis.
#'                             Default is 1
#'
#' @param x_coord_adjust Horizontal adjustment of the arrow origin, default is 0
#'
#' @param y_coord_adjust Vertical adjustment of the arrow origin, default is 0
#'
#' @param arrow_size The size of the arrow head, default is 0.2 cm
#'
#' @param arrow_width The width of the arrow, default is 0.9
#'
#' @param arrow_text_color The color of the arrow and label, default is "black"
#'
#' @param arrow_type The type of the arrow, default is "closed"
#'
#' @param text_font_size The font size of the text, default is 4
#'
#' @param relative_text_dst The relative distance of the text from the arrow,
#'                          default is 0.15
#'
#' @param ... Additional arguments to pass to the SCpubr::do_DimPlot function
#'
#' @return a ggplot object with the dimension arrows added
#'
#' @import SCpubr
#'
#' @export
#'
do_DimPlot_arrows <- function(
  object,
  reduction = "umap",
  dims = NULL,
  group.by = NULL,
  colors.use = NULL,
  pt.size = 0.5,
  shuffle = FALSE,
  label = FALSE,
  repel = label,
  legend.position = "right",
  arrow_length = 0.1,
  arrow_x_adjust_ratio = 1,
  arrow_y_adjust_ratio = 1,
  x_coord_adjust = 0,
  y_coord_adjust = 0,
  arrow_size = 0.2,
  arrow_width = 0.9,
  arrow_text_color = "black",
  arrow_type = "closed",
  text_font_size = 4,
  relative_text_dst = 0.15,
  ...
) {
  # Generate the plot
  plot <- SCpubr::do_DimPlot(
    object,
    reduction = reduction,
    group.by = group.by,
    colors.use = colors.use,
    pt.size = pt.size,
    shuffle = shuffle,
    label = label,
    repel = repel,
    legend.position = legend.position,
    ...
  )

  # Add dimension arrows
  plot <- add_dimension_arrows(
    plot,
    reduction = reduction,
    dims = dims,
    arrow_length = arrow_length,
    arrow_x_adjust_ratio = arrow_x_adjust_ratio,
    arrow_y_adjust_ratio = arrow_y_adjust_ratio,
    x_coord_adjust = x_coord_adjust,
    y_coord_adjust = y_coord_adjust,
    arrow_size = arrow_size,
    arrow_width = arrow_width,
    arrow_text_color = arrow_text_color,
    arrow_type = arrow_type,
    text_font_size = text_font_size,
    relative_text_dst = relative_text_dst
  )

  return(plot)
}



#' @title do_FeaturePlot with arrows on the bottom left corner
#' 
#' @description This function is a further implementation of the do_FeaturePlot
#'             function from the SCpubr package. It generates a ggplot object
#'            with the dimension arrows added to the bottom left corner of
#'            the plot. The function is useful for visualizing single cell
#'           data with dimension reduction methods such as UMAP, t-SNE, PCA,
#'          and spatial coordinates.
#' 
#' @param object A Seurat object
#' 
#' @param features A vector of feature names to plot
#' 
#' @param assay The assay to use for plotting, default is the active assay
#' 
#' @param reduction The reduction method used to generate the plot, default is
#'                 "umap".
#' 
#' @param slot The slot to use for the reduction, default is "data"
#' 
#' @param order Whether to order the cells based on feature expression. Default
#'              is TRUE
#' 
#' @param dims A vector of dimension names, overwrites the default dimension
#'             names assumed from reduction method. (Overwrite the dims argument
#'             in the SCpubr::do_FeaturePlot function)
#' 
#' @param pt.size The size of the points, default is 0.5
#' 
#' @param legend.position The position of the legend, default is 'right'
#' 
#' @param arrow_length The size of the arrow, default is 0.15 (15% of the plot)
#' 
#' @param arrow_x_adjust_ratio The x-axis adjustment ratio for the arrow,
#' 
#' @param arrow_y_adjust_ratio The y-axis adjustment ratio for the arrow,
#' 
#' @param x_coord_adjust Horizontal adjustment of the arrow origin, default is 0
#' 
#' @param y_coord_adjust Vertical adjustment of the arrow origin, default is 0
#' 
#' @param arrow_size The size of the arrow head, default is 0.2 cm
#' 
#' @param arrow_width The width of the arrow, default is 0.9
#' 
#' @param arrow_text_color The color of the arrow and label, default is "black"
#' 
#' @param arrow_type The type of the arrow, default is "closed"
#' 
#' @param text_font_size The font size of the text, default is 4
#' 
#' @param relative_text_dst The relative distance of the text from the arrow,
#' 
#' @param ... Additional arguments to pass to the SCpubr::do_FeaturePlot function
#' 
#' @return a ggplot object with the dimension arrows added
#' 
#' @import SCpubr
#' 
#' @export
#' 
do_FeaturePlot_arrows <- function(
  object,
  features,
  assay = NULL,
  reduction = "umap",
  slot = "data",
  order = TRUE,
  dims = NULL,
  pt.size = 0.5,
  legend.position = "right",
  arrow_length = 0.1,
  arrow_x_adjust_ratio = 1,
  arrow_y_adjust_ratio = 1,
  x_coord_adjust = 0,
  y_coord_adjust = 0,
  arrow_size = 0.2,
  arrow_width = 0.9,
  arrow_text_color = "black",
  arrow_type = "closed",
  text_font_size = 4,
  relative_text_dst = 0.15,
  ...
) {
  # Generate the plot
  plot <- SCpubr::do_FeaturePlot(
    object,
    features = features,
    assay = assay,
    reduction = reduction,
    slot = slot,
    order = order,
    dims = dims,
    pt.size = pt.size,
    legend.position = legend.position,
    ...
  )

  # Add dimension arrows
  plot <- add_dimension_arrows(
    plot,
    reduction = reduction,
    dims = dims,
    arrow_length = arrow_length,
    arrow_x_adjust_ratio = arrow_x_adjust_ratio,
    arrow_y_adjust_ratio = arrow_y_adjust_ratio,
    x_coord_adjust = x_coord_adjust,
    y_coord_adjust = y_coord_adjust,
    arrow_size = arrow_size,
    arrow_width = arrow_width,
    arrow_text_color = arrow_text_color,
    arrow_type = arrow_type,
    text_font_size = text_font_size,
    relative_text_dst = relative_text_dst
  )

  return(plot)
}