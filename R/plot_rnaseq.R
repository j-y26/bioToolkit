# Simplified Plotting for RNA-Seq Data

#' @title Volcano plot for gene expression analysis
#'
#' @description Generate a volcano plot for gene expression analysis based on
#'              the log2 fold change and the -log10 p-value.
#' 
#' @details This function inherently uses the EnhancedVolcano package to
#'          generate a volcano plot for gene expression analysis. Some
#'          aesthetic are modified from the default plot to make it more
#'          suitable for RNA-seq data analysis. Use the \code{...} to provide
#'          additional arguments to EnhanceVolcano::EnhancedVolcano.
#'
#' @param data A data frame containing the gene expression data
#'
#' @param de_method The method used for differential expression analysis.
#'                  Must be one of "edgeR" or "DESeq2". Default is "edgeR".
#'                  Use this argument to automatically set the column names
#'                  for log2 fold change and adjusted P-value. If provided value
#'                  is not "edgeR" or "DESeq2", the user must provide both the
#'                  column names for log2 fold change and P-value using the
#'                  \code{x} and \code{y} arguments, respectively.
#'
#' @param x The column name of the log2 fold change. If not \code{NULL},
#'          overrides the value obtained from de_method.
#'
#' @param y The column name of the P-value for use. If not \code{NULL},
#'          overrides the value obtained from de_method.
#'
#' @param label The column name of the gene names
#'
#' @param selected_labels A vector of gene names to be highlighted in the plot
#'
#' @param p_cutoff The p-value cutoff for significance
#'
#' @param fc_cutoff The log2 fold change cutoff for significance
#'
#' @param up_color The color for up-regulated genes
#'
#' @param down_color The color for down-regulated genes
#'
#' @param non_de_color The color for non-differentially expressed genes
#'
#' @param up_label The label for up-regulated genes
#'
#' @param down_label The label for down-regulated genes
#'
#' @param non_de_label The label for non-differentially expressed genes
#'
#' @param number_label A logical value indicating whether to display the number
#'                     of up and down-regulated genes
#' 
#' @param num_label_vjust The vertical justification for the number labels
#'
#' @param num_label_size The size of the number labels
#'
#' @param x_lim The limits for the x-axis
#'
#' @param y_lim The limits for the y-axis
#'
#' @param x_lab The label for the x-axis
#'
#' @param y_lab The label for the y-axis
#'
#' @param axis_lab_size The size of the axis labels
#'
#' @param lab_size The size of the gene labels
#'
#' @param bold_text A logical value indicating whether to bold the text
#'
#' @param legend_lab_size The size of the legend labels
#'
#' @param legend_lab_position The position of the legend labels
#'
#' @param legend_icon_size The size of the legend icons
#'
#' @param draw_connectors A logical value indicating whether to draw connectors
#'
#' @param max_overlaps The maximum number of overlaps for the gene labels
#'
#' @param ... Additional arguments to be passed to the plot function
#'
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @import ggplot2 dplyr
#'
#' @return A volcano plot
#'
#' @export
#'
plot_volcano <- function(
    data,
    de_method = "edgeR",
    x = NULL,
    y = NULL,
    label = rownames(data),
    selected_labels = NULL,
    p_cutoff = 0.05,
    fc_cutoff = 0,
    up_color = "firebrick4",
    down_color = "midnightblue",
    non_de_color = "grey",
    up_label = "Up-regulated",
    down_label = "Down-regulated",
    non_de_label = "Non-DE",
    number_label = TRUE,
    num_label_vjust = 0,
    num_label_size = 4.5,
    x_lim = c(
      min(data[[x]], na.rm = TRUE) - 0.5,
      max(data[[x]], na.rm = TRUE) + 0.5
    ),
    y_lim = c(0, max(-log10(data[[y]]), na.rm = TRUE) + 1),
    x_lab = bquote(~ log[2] ~ "fold change"),
    y_lab = bquote(~ -log[10] ~ adj.P),
    axis_lab_size = 12,
    lab_size = 4,
    bold_text = TRUE,
    legend_lab_size = 12,
    legend_lab_position = "top",
    legend_icon_size = 5,
    draw_connectors = TRUE,
    max_overlaps = 10,
    ...) {
  # Check if the DE method is valid, case-insensitive
  de_method <- tolower(de_method)
  if (de_method == "edger" && is.null(x) && is.null(y)) {
    x <- FC_COL_EDGER
    y <- PVAL_COL_EDGER
  } else if (de_method == "deseq2" && is.null(x) && is.null(y)) {
    x <- FC_COL_DESEQ2
    y <- PVAL_COL_DESEQ2
  } else {
    if (is.null(x) || is.null(y)) {
      stop("Please provide the column names for log2 fold change and P-value.")
    }
  }

  # Define key-value pairs for color coding
  key_val <- ifelse(
    data[[y]] < p_cutoff & data[[x]] > fc_cutoff, up_color,
    ifelse(
      data[[y]] < p_cutoff & data[[x]] < -fc_cutoff, down_color,
      non_de_color
    )
  )
  key_val[is.na(key_val)] <- non_de_color
  names(key_val)[key_val == up_color] <- up_label
  names(key_val)[key_val == down_color] <- down_label
  names(key_val)[key_val == non_de_color] <- non_de_label

  # Calculate the number of up and down-regulated genes
  num_up <- sum(data[[y]] < p_cutoff & data[[x]] > fc_cutoff) %>% as.character()
  num_down <- sum(data[[y]] < p_cutoff & data[[x]] < -fc_cutoff) %>% as.character()

  # Parse whether bold texts are used
  if (bold_text) {
    lab_face <- "bold"
    x_lab <- bquote(bold(.(x_lab)))
    y_lab <- bquote(bold(.(y_lab)))
    num_up <- bquote(bold(.(num_up)))
    num_down <- bquote(bold(.(num_down)))
  } else {
    lab_face <- "plain"
  }

  # Create the volcano plot
  v_plot <- EnhancedVolcano::EnhancedVolcano(
    data,
    lab = label,
    selectLab = selected_labels,
    labSize = lab_size,
    x = x,
    y = y,
    title = "",
    subtitle = "",
    caption = "",
    pCutoff = p_cutoff,
    FCcutoff = fc_cutoff,
    colCustom = key_val,
    xlim = x_lim,
    ylim = y_lim,
    xlab = x_lab,
    ylab = y_lab,
    axisLabSize = axis_lab_size,
    legendLabSize = legend_lab_size,
    legendPosition = legend_lab_position,
    legendIconSize = legend_icon_size,
    labFace = lab_face,
    drawConnectors = draw_connectors,
    max.overlaps = max_overlaps,
    ...
  )

  # Add the number of up and down-regulated genes at the bottom right and
  # left corners, respectively
  if (number_label) {
    v_plot <- v_plot + annotate(
      "text",
      x = x_lim[2],
      y = y_lim[1],
      label = num_up,
      hjust = 1,
      vjust = num_label_vjust,
      size = num_label_size
    ) + annotate(
      "text",
      x = x_lim[1],
      y = y_lim[1],
      label = num_down,
      hjust = 0,
      vjust = num_label_vjust,
      size = num_label_size
    )
  }

  # Bold the rest of the text if bold_text is TRUE
  if (bold_text) {
    v_plot <- v_plot + theme(
      axis.text = element_text(face = "bold"),
      axis.title.x = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      text = element_text(face = "bold")
    )
  }

  return(v_plot)
}
