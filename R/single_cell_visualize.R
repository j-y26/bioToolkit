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