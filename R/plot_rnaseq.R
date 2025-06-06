# Simplified Plotting for RNA-Seq Data

#' @title Volcano plot for gene expression analysis
#'
#' @description Generate a volcano plot for gene expression analysis based on
#'              the log2 fold change and the -log10 p-value.
#'
#' @details This function generates a volcano plot for gene expression analysis.
#'
#' @param data A data frame containing the gene expression data
#'
#' @param de_method The method used for differential expression analysis.
#'                  Must be one of "edgeR", "DESeq2", "limma", or "seurat". Default is "edgeR".
#'                  Use this argument to automatically set the column names
#'                  for log2 fold change and adjusted P-value. If provided value
#'                  is not "edgeR", "DESeq2", "limma", or "seurat", the user must provide both the
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
#' @param draw_connectors A logical value indicating whether to draw connectors
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
#' @param alpha The transparency level of the points in the plot
#'
#' @param pt_size The size of the points in the plot
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
#' @param legend_lab_size The size of the legend labels
#'
#' @param legend_lab_position The position of the legend labels
#'
#' @param legend_icon_size The size of the legend icons
#'
#' @param max_overlaps The maximum number of overlaps for the gene labels
#'
#' @import ggplot2 dplyr rlang
#'
#' @importFrom ggrepel geom_text_repel
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
    label = "row.names",
    selected_labels = rownames(data),
    draw_connectors = TRUE,
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
    alpha = 0.5,
    pt_size = 2,
    x_lim = NULL,
    y_lim = NULL,
    x_lab = bquote(~ log[2] ~ "fold change"),
    y_lab = bquote(~ -log[10] ~ adj.P),
    axis_lab_size = 12,
    lab_size = 3.5,
    legend_lab_size = 12,
    legend_lab_position = "top",
    legend_icon_size = 4,
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
  } else if (de_method == "limma" && is.null(x) && is.null(y)) {
    x <- FC_COL_LIMMA
    y <- PVAL_COL_LIMMA
  } else if (de_method == "seurat" && is.null(x) && is.null(y)) {
    x <- FC_COL_SEURAT
    y <- PVAL_COL_SEURAT
  } else {
    if (is.null(x) || is.null(y)) {
      stop("Please provide the column names for log2 fold change and P-value.")
    }
  }

  # Add additional information for plotting
  if (label == "row.names") {
    data <- data %>% mutate(label = rownames(data))
  } else {
    if (!label %in% colnames(data)) {
      stop("The label column does not exist in the data frame.")
    } else {
      data <- data %>% mutate(label = .data[[label]])
    }
  }

  data <- data %>%
    mutate(
      !!y := ifelse(.data[[y]] == 0, 1e-300, .data[[y]]),
      log_p = -log10(.data[[y]]),
      logfc = .data[[x]],
      group = case_when(
        .data[[y]] < p_cutoff & .data[[x]] > fc_cutoff ~ up_label,
        .data[[y]] < p_cutoff & .data[[x]] < -fc_cutoff ~ down_label,
        TRUE ~ non_de_label
      ),
      color = case_when(
        group == up_label ~ up_color,
        group == down_label ~ down_color,
        TRUE ~ non_de_color
      )
    ) %>%
      mutate(group = factor(group, levels = c(up_label, down_label, non_de_label)))

  # Calculate the number of up and down-regulated genes
  num_up <- sum(data[[y]] < p_cutoff & data[[x]] > fc_cutoff) %>% as.character()
  num_down <- sum(data[[y]] < p_cutoff & data[[x]] < -fc_cutoff) %>% as.character()

  # Define the limits for the x and y axes
  if (is.null(x_lim)) {
    x_lim <- c(min(data[[x]], na.rm = TRUE) - 0.5, max(data[[x]], na.rm = TRUE) + 0.5)
  }
  if (is.null(y_lim)) {
    y_lim <- c(0, max(data$log_p, na.rm = TRUE) + 1)
  }

  # Color assignment
  named_colors <- setNames(
    c(up_color, down_color, non_de_color),
    c(up_label, down_label, non_de_label)
  )

  # Create the volcano plot
  v_plot <- ggplot(data, aes(x = logfc, y = log_p)) +
    geom_point(aes(color = group), size = pt_size, alpha = alpha, stroke = NA) +
    scale_color_manual(values = named_colors, drop = TRUE) +
    labs(x = x_lab, y = y_lab) +
    xlim(x_lim) + ylim(y_lim) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black", linewidth = 0.3) +
    theme_classic(base_size = axis_lab_size) +
    theme(
      legend.position = legend_lab_position,
      legend.title = element_blank(),
      legend.text = element_text(size = legend_lab_size),
    )

  if (fc_cutoff > 0) {
      v_plot <- v_plot +
          geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
                     linetype = "dashed", color = "black", linewidth = 0.3)
  }

  # Prepare the labels for the selected genes
  if (!is.null(selected_labels)) {
    # Here we use having to label more than half or 100 genes, whichever is bigger
    # to avoid overplotting

    # If less than half or 100 genes, use the selected labels
    if (length(selected_labels) < max(100, nrow(data) / 2)) {
      label_data <- data %>%
        filter(label %in% selected_labels) %>%
        select(label, logfc, log_p)
    } else {
      # Otherwise, choose the top genes based on log_p and logfc
      # Use a score that is the product of log_p and logfc
      # Rank and choose the top 10 genes in each direction
      label_data <- data %>%
        filter(group %in% c(up_label, down_label)) %>%
        mutate(score = log_p * abs(logfc)) %>%
        arrange(desc(score)) %>%
        group_by(group) %>%
        slice_head(n = 10) %>%
        ungroup() %>%
        select(label, logfc, log_p)
    }
  } else {
    label_data <- NULL
  }

  # Add the labels to the plot
  if (!is.null(label_data)) {
    if (draw_connectors) {
      v_plot <- v_plot + ggrepel::geom_text_repel(
        data = label_data,
        aes(x = logfc, y = log_p, label = label),
        size = lab_size,
        max.overlaps = max_overlaps,
        segment.color = "black",
        segment.size = 0.5,
        segment.alpha = 0.5,
        nudge_x = 0.2,
        nudge_y = 0.2,
        fontface = "bold",
      )
    } else {
      v_plot <- v_plot + geom_text(
        data = label_data,
        aes(x = logfc, y = log_p, label = label),
        size = lab_size,
        nudge_x = 0.2,
        nudge_y = 0.2,
        fontface = "bold"
      )
    }
  }


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
      size = num_label_size,
      fontface = "bold"
    ) + annotate(
      "text",
      x = x_lim[1],
      y = y_lim[1],
      label = num_down,
      hjust = 0,
      vjust = num_label_vjust,
      size = num_label_size,
      fontface = "bold"
    )
  }

  # Increase the size of the legend icons
  v_plot <- v_plot +
      guides(color = guide_legend(override.aes = list(size = legend_icon_size)))

  return(v_plot)
}


#' @title Plotting top PCA loadings
#'
#' @description This function extracts the PCA loadings from a PCA object
#'              generated using the \code{prcomp} function and plots the top
#'              genes for each principal component.
#'
#' @param pca_obj A PCA object generated using the \code{prcomp} function.
#'
#' @param pcs A vector of integers specifying the principal components to
#'            extract the loadings. Default is to extract PC 1-3.
#'
#' @param top_n The number of top genes to extract for each principal component.
#'              Default is 20.
#'
#' @param low_color The color for low loadings
#'
#' @param mid_color The color for mid loadings
#'
#' @param high_color The color for high loadings
#'
#' @export
#'
#' @return A list, where each element is a PC, and each element is a named vector
#'         containing the top genes and their loadings.
#'
#' @importFrom ggplot2 ggplot geom_bar coord_flip facet_wrap scale_fill_gradient2
#' @importFrom dplyr arrange group_by ungroup
#'
plot_top_pca_loadings <- function(pca_obj,
                                  pcs = 1:3,
                                  top_n = 20,
                                  low_color = "darkgreen",
                                  mid_color = "white",
                                  high_color = "darkslateblue") {
  top_loadings <- extract_pc_loadings(pca_obj, pcs, top_n)

  # Convert list to data frame
  plot_data <- do.call(rbind, lapply(names(top_loadings), function(pc) {
    data.frame(Gene = names(top_loadings[[pc]]),
               Loading = unname(top_loadings[[pc]]),
               PC = pc)
  }))

  # Arrange from most positive to most negative
  plot_data <- plot_data %>%
    mutate(PC = gsub("_", "-", PC)) %>%
    group_by(PC) %>%
    arrange(desc(Loading), .by_group = TRUE) %>%
    ungroup()

  # Plot
  p <- ggplot(plot_data, aes(x = reorder(Gene, Loading), y = Loading, fill = Loading)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~ PC, scales = "free_y") +
    scale_fill_gradient2(low = low_color, mid = mid_color, high = high_color, midpoint = 0) +
    theme_minimal() +
    labs(title = "Top Contributing Genes to Principal Components",
         x = "Gene",
         y = "Loading",
         fill = "Direction") +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          strip.text = element_text(size = 12, face = "bold"))

  return(p)
}
