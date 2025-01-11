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
#' @param legend.framewidth The width of the legend frame, default is 0.3
#' 
#' @param legend.tickwidth The width of the legend ticks, default is 0.3
#' 
#' @param legend.length The length of the legend, default is 10
#' 
#' @param legend.width The width of the legend, default is 0.6
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
#' @param sequential.palette The color palette to use for the feature plot
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
  sequential.palette = "Purples",
  legend.position = "right",
  legend.framewidth = 0.3,
  legend.tickwidth = 0.3,
  legend.length = 10,
  legend.width = 0.6,
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
    legend.framewidth = legend.framewidth,
    legend.tickwidth = legend.tickwidth,
    legend.length = legend.length,
    legend.width = legend.width,
    sequential.palette = sequential.palette,
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
#' @return a ggplot object with the dot plot
#' 
#' @import ggplot2 dplyr
#' @importFrom RColorBrewer brewer.pal
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
  sig_only = FALSE,
  invert = FALSE
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

  # Create the plot
  if (!invert) {
    dot_plot <- ggplot(de_results, aes(y = !!sym(comparison_var), x = !!sym(feature_var))) +
      geom_point(aes(size = dot_size, color = ifelse(signif, !!sym(logFC_var), NA)))
  } else {
    dot_plot <- ggplot(de_results, aes(y = !!sym(feature_var), x = !!sym(comparison_var))) +
      geom_point(aes(size = dot_size, color = ifelse(signif, !!sym(logFC_var), NA)))
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
      name = "log2(Fold Change)"
    ) +
    scale_size_continuous(name = "-log10(p-value)") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )

  return(dot_plot)
}