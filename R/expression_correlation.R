# Explore the correlation between gene expression and another variable

#' @title Calculate the correlation between gene expression and another variable
#'        using a linear model
#'
#' @description This function calculates the correlation between gene expression
#'              and another continuous variable, most commonly another gene.
#'              This function is used to identify the list of genes in which the
#'              expression is correlated with another variable/gene in a dose-
#'              dependent manner. This function uses a linear model to calculate
#'              the correlation, with the query gene as the independent variable
#'              and all other genes as the dependent variable. The significance
#'              of the correlation is calculated using the p-value of the linear
#'              model and is adjusted for multiple testing using a user-defined
#'              method. The slope of the linear model is used as the correlation
#'              coefficient, with R-squared as the measure of goodness of fit.
#'              The strength of the correlation is calculated using a
#'              user-defined method (e.g., Pearson, Spearman, Kendall).
#'
#' @param gene_expr A data frame containing the gene expression data. The data
#'                  frame should have the genes as rows and samples as columns.
#'                  The gene expression data should be normalized (e.g. CPM)
#'                  and consist of only numeric values. All genes in this data
#'                  frame, except the query gene, will be used as the dependent
#'                  variable in the linear model. However, the query_var can
#'                  also be specified in the \code{metadata} data frame. Note,
#'                  only numeric values are allowed in the gene expression data.
#'
#' @param query_var A character string specifying the query gene for which the
#'                  correlation with other genes will be calculated. The query
#'                  gene should be present in the gene expression data frame,
#'                  and it must be a row name in the \code{gene_expr} data
#'                  frame. In general, the query gene is the gene of interest
#'                  and should be a continuous variable specified in the
#'                  \code{gene_expr} data frame. However, in other cases, the
#'                  query gene can be a categorical/ordinal variable specified
#'                  in the \code{metadata} data frame. In this case, any
#'                  ordinal variable should be converted into a numeric variable
#'                  before using this function. Categorical variables should be
#'                  converted into a factor and have the base level set as the
#'                  reference level.
#'
#' @param metadata A data frame containing the metadata information for the
#'                 samples. The data frame should have samples as rows and
#'                 metadata variables as columns. The metadata information can
#'                 contain the \code{query_var} and/or \code{vars_to_regress}
#'                 variables that are used as covariates in the linear model.
#'                 The metadata information can be used to adjust the gene
#'                 expression data for confounding variables. The metadata
#'                 can contain both categorical and continuous variables.
#'                 ordinal variables should be converted into numeric variables
#'                 before using this function. Categorical variables should be
#'                 converted into a factor and have the base level set as the
#'                 reference level. Default value is \code{NULL}. Note, row
#'                 names of the \code{metadata} data frame should match the
#'                 column names of the \code{gene_expr} data frame.
#'
#' @param vars_to_regress A character vector specifying the variables to regress
#'                        out from the linear model. The default value is
#'                        \code{NULL}, which means no variables will be
#'                        regressed out. The variables can present in either
#'                        the \code{gene_expr} or \code{metadata} data frame.
#'                        However, if the variable are not continuous, they
#'                        should only exist in the \code{metadata} data frame.
#'
#' @param scale_expr A logical value indicating whether to scale the gene
#'                   expression data before calculating the correlation. If
#'                   \code{TRUE}, the gene expression data will be scaled using
#'                   the \code{scale} function to a mean of 0 and a standard
#'                   deviation of 1 (z-score). Expression of each gene will be
#'                   scaled independently across samples. The default value is
#'                   \code{TRUE}. If false, the original gene expression data
#'                   will be used for calculations.
#'
#' @param plot_model A logical value indicating whether to plot the linear model
#'                   for each gene. The default value is \code{FALSE}. If true,
#'                   a plot will be generated for each gene showing the linear
#'                   model with the query gene on the x-axis and the gene
#'                   expression on the y-axis. The plot will also show the
#'                   linear model coefficients, R-squared value, correlation
#'                   coefficient, and FDR-adjusted p-value. The plot will be
#'                   a ggplot object and can be further customized using the
#'                   ggplot2 functions. The plot will be saved in the
#'                   \code{model_plots_positive} and \code{model_plots_negative}
#'                   directories in the \code{base_dir} directory.
#'
#' @param base_dir A character string specifying the base directory to save the
#'                 model results and plots. The default value is the current
#'                 working directory. The model results will be saved as a CSV
#'                 file named \code{query_var_correlation_results.csv}.
#'
#' @param corr_method A character string specifying the method to calculate the
#'                    correlation. The available methods include "pearson",
#'                    "spearman", "kendall". The default method is "pearson".
#'
#' @param p_adjust_method A character string specifying the method to adjust the
#'                        p-values for multiple testing. The available methods
#'                        include "BH" (Benjamini-Hochberg), "bonferroni",
#'                        "holm", "hochberg", "hommel", "BY", "fdr", "none".
#'                        The default method is "BH". See \code{\link{p.adjust}}
#'                        for more details.
#'
#' @param file_type A character string specifying the file type for the model
#'                  plots. The available file types include "png", "jpeg",
#'                  "tiff", "pdf", "svg". The default file type is "png".
#'
#' @param dpi An integer specifying the resolution of the model plots in dots
#'            per inch (DPI). The default value is 300. The higher the DPI,
#'            the better the quality of the plots. The DPI value should be
#'            greater than 0.
#'
#' @param plot_group_var A character string specifying the variable in the
#'                       metadata data frame to group the samples in the plot.
#' 
#' @param switch_axes A logical value indicating whether to switch the axes of
#'                    the plot. The default value is \code{FALSE}. If
#'                    \code{TRUE}, the query gene will be plotted on the y-axis
#'                    and the gene expression on the x-axis.
#'
#' @param alpha The significance level for the adjusted p-values. The default
#'              value is 0.05.
#'
#' @param verbose A logical value indicating whether to print the progress
#'                messages. The default value is TRUE.
#'
#' @import dplyr tibble
#'
#' @export
#'
#' @return A list of four elements:
#'        \item{model_list}{A list of linear models for each gene, with the
#'                          query gene as the independent variable and the gene
#'                          expression as the dependent variable.}
#'        \item{model_results}{A data frame containing the correlation results,
#'                             including the correlation coefficient, p-value,
#'                             adjusted p-value, and R-squared value for each
#'                             gene.}
#'        \item{p_volcano}{A ggplot object of the volcano plot showing the
#'                         correlation results.}
#'        \item{pos_plots}{A list containing the ggplot objects for each gene
#'                           showing the linear model plot. The plots will be
#'                           saved automatically. Only available if
#'                          \code{plot_model} is \code{TRUE}.}
#'        \item{neg_plots}{A list containing the ggplot objects for each gene
#'                          showing the linear model plot. The plots will be
#'                          saved automatically. Only available if
#'                          \code{plot_model} is \code{TRUE}.}
#'
calc_expr_corr <- function(
  gene_expr,
  query_var,
  metadata = NULL,
  vars_to_regress = NULL,
  scale_expr = TRUE,
  plot_model = FALSE,
  base_dir = getwd(),
  corr_method = "pearson",
  p_adjust_method = "BH",
  file_type = "png",
  dpi = 300,
  plot_group_var = NULL,
  switch_axes = FALSE,
  alpha = 0.05,
  verbose = TRUE
) {


  # Check if the correlation method is valid
  corr_method <- tolower(corr_method)
  if (!corr_method %in% c("pearson", "spearman", "kendall")) {
    stop("Invalid correlation method. Please choose from 'pearson', 'spearman', 'kendall'.")
  }

  # Check if the p-value adjustment method is valid
  if (!p_adjust_method %in% c("BH", "bonferroni", "holm", "hochberg", "hommel", "BY", "fdr", "none")) {
    stop("Invalid p-value adjustment method. See ?p.adjust for available methods.")
  }

  # Check if plots too be plotted
  if (plot_model) {
    if (!dir.exists(file.path(base_dir))) {
      dir.create(base_dir)
    }
    if (!dir.exists(file.path(base_dir, "model_plots_positive"))) {
      dir.create(file.path(base_dir, "model_plots_positive"))
    }
    if (!dir.exists(file.path(base_dir, "model_plots_negative"))) {
      dir.create(file.path(base_dir, "model_plots_negative"))
    }

    # Check file type and DPI
    if (!file_type %in% c("png", "jpeg", "tiff", "pdf", "svg")) {
      stop("Invalid file type. Please choose from 'png', 'jpeg', 'tiff', 'pdf', 'svg'.")
    }
    if (dpi <= 0) {
      stop("DPI should be greater than 0.")
    }
  }

  # If metadata is provided, check if it is a data frame and has row names
  # matching the column names of gene expression data
  if (!is.null(metadata)) {
    if (!is.data.frame(metadata)) {
      stop("Metadata should be a data frame.")
    }
    if (!all(colnames(gene_expr) %in% rownames(metadata))) {
      stop("Row names of metadata should match the column names of gene expression data.")
    }
    # Order the metadata rows to match the gene expression columns
    metadata <- metadata[colnames(gene_expr), ]
  }

  # Check if the variables are in the data
  if (!is.null(metadata)) {
    all_vars <- c(rownames(gene_expr), colnames(metadata))
  } else {
    all_vars <- rownames(gene_expr)
  }
  if (!query_var %in% all_vars) {
    stop("Query gene not found in the data provided.")
  }
  if (!is.null(vars_to_regress)) {
    if (!all(vars_to_regress %in% all_vars)) {
      stop("Variables to regress out not found in the data provided.")
    }
  }
  if (plot_model && !is.null(plot_group_var)) {
    if (!plot_group_var %in% colnames(metadata)) {
      stop("Group variable not found in the metadata data frame.")
    }
  }

  # Scale the gene expression data if required
  if (scale_expr) {
    gene_expr <- t(scale(t(gene_expr)))
    if (verbose) {
      cat("Gene expression data scaled.\n")
    }
  } else {
    gene_expr <- as.matrix(gene_expr)
    if (verbose) {
      cat("Using original gene expression data.\n")
    }
  }

  # Check if the gene expression data is numeric
  if (!is.numeric(gene_expr)) {
    stop("Gene expression data should be numeric.")
  }

  # Determine where the query gene is located
  if (query_var %in% rownames(gene_expr)) {
    query_loc <- "expression"
  } else {
    query_loc <- "metadata"
  }

  # Genes to fit
  genes_to_fit <- rownames(gene_expr) %>% setdiff(query_var)

  # Initialize a list of models
  model_list <- list()
  model_results <- list()

  # Create a data frame with the variables to regress out
  regress_data <- list()
  if (!is.null(vars_to_regress)) {
    for (var in vars_to_regress) {
      if (var %in% rownames(gene_expr)) {
        regress_data[[var]] <- gene_expr[var, , drop = FALSE]
      } else {
        regress_data[[var]] <- metadata[, var, drop = FALSE]
      }
    }
  }
  if (length(regress_data) > 0) {
    regress_data <- do.call(cbind, regress_data)
  }

  # Loop through each gene and fit a linear model
  if (verbose) {
    cat(paste0("Fitting linear models and calculating correlation using ", corr_method, " method.\n"))
    if (!is.null(vars_to_regress)) {
      cat("Regressing out variables: ", paste(vars_to_regress, collapse = ", "), "\n")
    }
  }

  # Total number of genes to fit
  total_genes <- length(genes_to_fit)
  pct_25 <- round(total_genes * 0.25)
  pct_50 <- round(total_genes * 0.50)
  pct_75 <- round(total_genes * 0.75)

  for (i in seq_len(total_genes)) {

    gene <- genes_to_fit[i]

    # Create a data frame with the gene expression and query gene
    if (query_loc == "expression") {
      query <- gene_expr[query_var, ]
    } else {
      query <- metadata[, query_var]
    }

    data <- data.frame(
      query = query,
      gene = gene_expr[gene, ]
    )

    if (!is.null(vars_to_regress)) {
      data <- cbind(data, regress_data)
    }

    # Fit a linear model
    if (!is.null(vars_to_regress)) {
      model <- lm(gene ~ query + ., data = data)
    } else {
      model <- lm(gene ~ query, data = data)
    }

    # Perform correlation test
    corr <- cor.test(data$query, data$gene, method = corr_method)

    # Store the model
    model_list[[gene]] <- model
    model_results[[gene]] <- list(
      gene = gene,
      query = query_var,
      p_value = summary(model)$coefficients[2, 4],
      coef = coef(model)[2],
      r_squared = summary(model)$r.squared,
      corr = corr$estimate,
      corr_p_value = corr$p.value,
      corr_method = corr_method
    )


    # Print progress messages
    if (verbose) {
      if (i == pct_25) {
        cat("25% of genes fitted.\n")
      } else if (i == pct_50) {
        cat("50% of genes fitted.\n")
      } else if (i == pct_75) {
        cat("75% of genes fitted.\n")
      } else if (i == total_genes) {
        cat("All gene models fitted.\n")
      }
    }
  }

  # Convert the model results to a data frame
  model_results <- model_results %>%
    map_dfr(~ .x) %>%
    as_tibble()

  # Adjust the p-values
  model_results$FDR <- p.adjust(model_results$p_value, method = p_adjust_method)
  model_results$corr_FDR <- p.adjust(model_results$corr_p_value, method = p_adjust_method)

  model_results <- model_results %>%
    select(gene, query, coef, r_squared, p_value, FDR, corr, corr_p_value, corr_FDR, corr_method) %>%
    arrange(FDR) %>%
    mutate(sig = ifelse(FDR < alpha & coef > 0, "Positively correlated",
                        ifelse(FDR < alpha & coef < 0, "Negatively correlated", "Not correlated")))

  write.csv(model_results, file.path(base_dir, paste0(query_var, "_correlation_results.csv")), row.names = FALSE)

  # Plot a volcano plot
  p_volcano <- plot_volcano(
    data = model_results,
    x = "corr",
    y = "FDR",
    label = model_results$gene,
    x_lab = paste0(corr_method, " correlation"),
    up_label = "Pos. corr.",
    down_label = "Neg. corr.",
    non_de_label = "Not corr."
  )

  # Save the volcano plot
  ggsave(file.path(base_dir, paste0(query_var, "_correlation_volcano.", file_type)),
         plot = p_volcano, dpi = dpi, width = 5, height = 6)

  results <- list(
    model_list = model_list,
    model_results = model_results,
    p_volcano = p_volcano
  )

  # Plot the linear model for each gene
  if (plot_model) {
    pos_plots <- list()
    neg_plots <- list()

    cat("Plotting linear models for each gene.\n")

    for (gene in model_results$gene) {
      if (model_results$sig[model_results$gene == gene] == "Not correlated") {
        next
      }

      # Get the model
      model <- model_list[[gene]]

      min_query <- min(model$model$query)
      max_query <- max(model$model$query)

      min_gene <- min(model$model$gene)
      max_gene <- max(model$model$gene)

      # Statistics
      FDR <- signif(model_results$FDR[model_results$gene == gene], 4)
      slope <- signif(coef(model)["query"], 4)
      r_squared <- signif(summary(model)$r.squared, 4)
      corr <- signif(model_results$corr[model_results$gene == gene], 4)

      if (switch_axes) {
        # Plot the gene expression against the query gene
        if (!is.null(plot_group_var)) {
          sample_groups <- metadata[rownames(model$model), plot_group_var]
          p_gene_query <- ggplot(model$model,
                                aes(x = gene, y = query, color = sample_groups)) +
            geom_point() +
            stat_smooth(method = "lm", col = "firebrick4") +
            labs(
              x = gene,
              y = query_var,
              color = element_blank()
            ) +
            theme_classic() +
            annotate("text", x = max_gene, y = max_query, label = bquote("FDR =" ~ .(FDR) ~ "\n" ~
                              "Slope =" ~ .(slope) ~ "\n" ~
                              "R"^2 ~ " =" ~ .(r_squared) ~ "\n" ~
                              "Corr =" ~ .(corr)), hjust = 1, vjust = 1)
        } else {
          sample_groups <- NULL
          p_gene_query <- ggplot(model$model,
                                aes(x = gene, y = query)) +
            geom_point() +
            stat_smooth(method = "lm", col = "firebrick4") +
            labs(
              x = gene,
              y = query_var
            ) +
            theme_classic() +
            annotate("text", x = max_gene, y = max_query, label = bquote("FDR =" ~ .(FDR) ~ "\n" ~
                              "Slope =" ~ .(slope) ~ "\n" ~
                              "R"^2 ~ " =" ~ .(r_squared) ~ "\n" ~
                              "Corr =" ~ .(corr)), hjust = 1, vjust = 1)
        }

      } else {
        # Plot the gene expression against the query gene
        if (!is.null(plot_group_var)) {
          sample_groups <- metadata[rownames(model$model), plot_group_var]
          p_gene_query <- ggplot(model$model,
                                aes(x = query, y = gene, color = sample_groups)) +
            geom_point() +
            stat_smooth(method = "lm", col = "firebrick4") +
            labs(
              x = query_var,
              y = gene,
              color = element_blank()
            ) +
            theme_classic() +
            annotate("text", x = max_query, y = max_gene, label = bquote("FDR =" ~ .(FDR) ~ "\n" ~
                              "Slope =" ~ .(slope) ~ "\n" ~
                              "R"^2 ~ " =" ~ .(r_squared) ~ "\n" ~
                              "Corr =" ~ .(corr)), hjust = 1, vjust = 1)
        } else {
          sample_groups <- NULL
          p_gene_query <- ggplot(model$model,
                                aes(x = query, y = gene)) +
            geom_point() +
            stat_smooth(method = "lm", col = "firebrick4") +
            labs(
              x = query_var,
              y = gene
            ) +
            theme_classic() +
            annotate("text", x = max_query, y = max_gene, label = bquote("FDR =" ~ .(FDR) ~ "\n" ~
                              "Slope =" ~ .(slope) ~ "\n" ~
                              "R"^2 ~ " =" ~ .(r_squared) ~ "\n" ~
                              "Corr =" ~ .(corr)), hjust = 1, vjust = 1)
        }
      }

      # Save the plot
      if (model_results$sig[model_results$gene == gene] == "Positively correlated") {
        pos_plots[[gene]] <- p_gene_query
        ggsave(file.path(base_dir, "model_plots_positive", paste0(gene, "_pos_", query_var, "_correlation.", file_type)),
               plot = p_gene_query, dpi = dpi, width = 6, height = 4)
      } else {
        neg_plots[[gene]] <- p_gene_query
        ggsave(file.path(base_dir, "model_plots_negative", paste0(gene, "_neg_", query_var, "_correlation.", file_type)),
               plot = p_gene_query, dpi = dpi, width = 6, height = 4)
      }
    }

    results$pos_plots <- pos_plots
    results$neg_plots <- neg_plots

    cat("Linear models plotted.\n")
  }

  return(results)
}
