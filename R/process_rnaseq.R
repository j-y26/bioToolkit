# Common functionalities for processing RNA-seq data

#' @title Calculate ranks of genes for differential expression analysis
#'
#' @description This function calculates ranks of genes based on the
#'              differential expression analysis results for RNA-seq data.
#'              The pi-value is calculated using the following formula:
#'              \deqn{pi = -log10(p) * sign(log2fc)}
#'              where \eqn{p} is the (adjusted) p-value and \eqn{log2fc} is the
#'              log2 fold change. The rank is then calculated based on the
#'              pi-value.
#'
#' @param de_results A data frame containing the results of differential
#'                   expression analysis. The data frame should have at least
#'                   a column for log2 fold change and a column for (adjusted)
#'                   p-value.
#'
#' @param de_method The method used for differential expression analysis.
#'                  Must be one of "edgeR" or "DESeq2". Default is "edgeR".
#'                  Use this argument to automatically set the column names
#'                  for log2 fold change and adjusted P-value. If provided value
#'                  is not "edgeR" or "DESeq2", the user must provide both the
#'                  column names for log2 fold change and P-value using the
#'                  \code{log2fc_col} and \code{pvalue_col} arguments,
#'                  respectively.
#'
#' @param pvalue_col The name of the column in \code{de_results} that contains
#'                   the p-values. If not NULL, overrides the default column
#'                   name based on \code{de_method}.
#'
#' @param log2fc_col The name of the column in \code{de_results} that contains
#'                   the log2 fold changes.
#'
#' @export
#'
#' @return A data frame with the addition of the following columns:
#'         \itemize{
#'          \item{piValue}{The pi-value of the gene}
#'          \item{rank}{The rank of the gene based on the pi-value}
#'         }
#'
rank_de_genes <- function(
    de_results,
    de_method = "edgeR",
    pvalue_col = NULL,
    log2fc_col = NULL) {
  # Set column names based on the method used, case-insensitive
  de_method <- tolower(de_method)
  if (de_method == "edger" && is.null(pvalue_col) && is.null(log2fc_col)) {
    pvalue_col <- PVAL_COL_EDGER
    log2fc_col <- FC_COL_EDGER
  } else if (de_method == "deseq2" && is.null(pvalue_col) && is.null(log2fc_col)) {
    pvalue_col <- PVAL_COL_DESEQ2
    log2fc_col <- FC_COL_DESEQ2
  } else if (is.null(pvalue_col) || is.null(log2fc_col)) {
    stop("Please provide both the column names for log2 fold change and p-value.")
  }

  # Calculate pi-value
  de_results$piValue <- -log10(de_results[[pvalue_col]]) * sign(de_results[[log2fc_col]])

  # Rank genes based on pi-value, from max to min
  de_results$rank <- rank(-de_results$piValue, ties.method = "first")

  # re-order the rows based on rank
  de_results <- de_results[order(de_results$rank), ]

  return(de_results)
}


#' @title Filter genes based on differential expression analysis results
#' 
#' @description This function filters genes based on the results of differential
#'              expression analysis for RNA-seq data. The user can specify the
#'              minimum log2 fold change and maximum adjusted p-value to filter
#'              the genes.
#' 
#' @param de_results A data frame containing the results of differential
#'                   expression analysis. The data frame should have at least
#'                   a column for log2 fold change and a column for (adjusted)
#'                   p-value.
#' 
#' @param log2fc_min The minimum log2 fold change to filter genes. Default is
#'                   0, which means no filtering based on log2 fold change.
#' 
#' @param pvalue_max The maximum adjusted p-value to filter genes. Default is
#'                   0.05.
#' 
#' @param de_method The method used for differential expression analysis.
#'                  Must be one of "edgeR" or "DESeq2". Default is "edgeR".
#'                  Use this argument to automatically set the column names
#'                  for log2 fold change and adjusted P-value. If provided value
#'                  is not "edgeR" or "DESeq2", the user must provide both the
#'                  column names for log2 fold change and P-value using the
#'                  \code{log2fc_col} and \code{pvalue_col} arguments,
#'                  respectively.
#' 
#' @param pvalue_col The name of the column in \code{de_results} that contains
#'                   the p-values. If not NULL, overrides the default column
#'                   name based on \code{de_method}.
#'                  
#' @param log2fc_col The name of the column in \code{de_results} that contains
#'                   the log2 fold changes.
#' 
#' @export
#' 
#' @return A list of two data frames:
#'        \itemize{
#'          \item{up_reg}{Data frame containing up-regulated genes}
#'          \item{down_reg}{Data frame containing down-regulated genes}
#'        }
#' 
filter_de_genes <- function(
    de_results,
    log2fc_min = 0,
    pvalue_max = 0.05,
    de_method = "edgeR",
    pvalue_col = NULL,
    log2fc_col = NULL) {
  # Set column names based on the method used, case-insensitive
  de_method <- tolower(de_method)
  if (de_method == "edger" && is.null(pvalue_col) && is.null(log2fc_col)) {
    pvalue_col <- PVAL_COL_EDGER
    log2fc_col <- FC_COL_EDGER
  } else if (de_method == "deseq2" && is.null(pvalue_col) && is.null(log2fc_col)) {
    pvalue_col <- PVAL_COL_DESEQ2
    log2fc_col <- FC_COL_DESEQ2
  } else if (is.null(pvalue_col) || is.null(log2fc_col)) {
    stop("Please provide both the column names for log2 fold change and p-value.")
  }

  # Filter genes based on log2 fold change and p-value
  up_reg <- de_results[de_results[[log2fc_col]] > log2fc_min &
                          de_results[[pvalue_col]] < pvalue_max, ]
  down_reg <- de_results[de_results[[log2fc_col]] < -log2fc_min &
                            de_results[[pvalue_col]] < pvalue_max, ]

  return(list(up_reg = up_reg, down_reg = down_reg))
}


#' @title Filter out genes with low expression based on minimum group size
#' 
#' @description This function filters out genes with low expression based on
#'              the minimum group size. The user can specify the minimum number
#'              of samples in which the gene should be expressed. Essentially,
#'              a gene is kept only if it has at least 1 count per million
#'              (CPM) in at least \code{min_group_size} samples.
#' 
#' @param raw_counts A data frame containing the raw counts of genes. The data
#'                   frame should have samples as columns and genes as rows.
#'                   All values should be non-negative integers.
#' 
#' @param min_group_size The minimum number of samples in which the gene should
#'                       be expressed. Default is 4.
#' 
#' @importFrom edgeR cpm DGEList
#' 
#' @export
#' 
#' @return A data frame containing the filtered raw counts
#' 
filter_low_expression_genes <- function(raw_counts, min_group_size = 4) {
  # Convert raw counts to DGEList object
  dge <- DGEList(counts = raw_counts)

  # Calculate counts per million (CPM)
  cpm_values <- cpm(dge)

  # Filter genes based on minimum group size
  keep_genes <- rowSums(cpm_values > 1) >= min_group_size
  filtered_counts <- raw_counts[keep_genes, ]

  return(filtered_counts)
}


#' @title Summarize group expression values for each gene
#' 
#' @description This function summarizes the expression values of genes for each
#'              group. The user can specify the method to summarize the values
#'              for each group. The default method is to calculate the mean of
#'              the expression values.
#' 
#' @param counts A data frame containing the expression values of genes. The
#'               data frame should have samples as columns and genes as rows.
#'               All values should be numeric.
#' 
#' @param group_info A data frame containing the group information of samples.
#'                  The data frame should have samples as rows and groups as
#'                  columns. One column should contain the information to group
#'                  the samples.
#' 
#' @param group_by The name of the column in \code{group_info} that contains the
#'                 information to group the samples. Row names must be the same
#'                 as the column names in \code{counts}.
#' 
#' @param method The method to summarize the expression values for each group.
#'               Must be one of "mean", "median", "min", "max", or "sum".
#'               Default is "mean".
#' 
#' @param pre_scale A logical value indicating whether to pre-scale the counts
#'                  before summarizing. Default is FALSE. Scaling will be done
#'                  using the \code{scale} function. Scaling is done across
#'                  samples for each gene.
#' 
#' @param post_scale A logical value indicating whether to post-scale the
#'                   summarized values. Default is FALSE. Scaling will be done
#'                   using the \code{scale} function. Scaling is done across
#'                   groups for each gene.
#' 
#' @export
#' 
#' @return A data frame containing the summarized expression values for each
#'         gene across samples in each group.
#' 
summarize_group_expression <- function(counts,
                                       group_info,
                                       group_by,
                                       method = "mean",
                                       pre_scale = FALSE,
                                       post_scale = FALSE) {

  # Check if method is valid
  valid_methods <- c("mean", "median", "min", "max", "sum")
  if (!method %in% valid_methods) {
    stop("Invalid method. Please choose one of: mean, median, min, max, sum.")
  }

  # Check sample names match
  if (!all(colnames(counts) %in% rownames(group_info))) {
    stop("Some sample names in counts are missing from group_info.")
  }

  # Check if counts is numeric data frame
  if (!is.data.frame(counts) || !all(sapply(counts, is.numeric))) {
    stop("counts must be a numeric data frame.")
  }

  # Get unique groups and samples in each group
  groups <- unique(group_info[[group_by]])
  grouped_samples <- split(rownames(group_info), group_info[[group_by]])

  # Initialize a data frame to store the summarized expression values
  summarized_values <- data.frame(matrix(NA, nrow = nrow(counts), ncol = length(groups)))
  colnames(summarized_values) <- groups
  rownames(summarized_values) <- rownames(counts)

  # Summarize expression values for each group
  for (group in groups) {
    samples <- intersect(grouped_samples[[group]], colnames(counts))  # Ensure valid samples
    if (length(samples) == 0) next  # Skip if no valid samples

    group_counts <- counts[, samples, drop = FALSE]  # Prevent vector conversion

    # Pre-scale counts if required
    if (pre_scale) {
      group_counts <- group_counts %>%
        t() %>%
        scale() %>%
        t() %>%
        as.data.frame()
    }

    if (method == "mean") {
      summarized_values[, group] <- rowMeans(group_counts, na.rm = TRUE)
    } else if (method == "median") {
      summarized_values[, group] <- apply(group_counts, 1, median, na.rm = TRUE)
    } else if (method == "min") {
      summarized_values[, group] <- apply(group_counts, 1, min, na.rm = TRUE)
    } else if (method == "max") {
      summarized_values[, group] <- apply(group_counts, 1, max, na.rm = TRUE)
    } else if (method == "sum") {
      summarized_values[, group] <- rowSums(group_counts, na.rm = TRUE)
    }
  }

  # Post-scale summarized values if required
  if (post_scale) {
    summarized_values <- summarized_values %>%
      t() %>%
      scale() %>%
      t() %>%
      as.data.frame()
  }

  return(summarized_values)
}

#' @title Extract PCA loadings from PCA object
#' 
#' @description This function extracts the PCA loadings from a PCA object
#'              generated using the \code{prcomp} function.
#' 
#' @param pca_obj A PCA object generated using the \code{prcomp} function.
#' 
#' @param pcs A vector of integers specifying the principal components to
#'            extract the loadings. Default is to extract PC 1-3.
#' 
#' @param top_n The number of top genes to extract for each principal component.
#'              Default is 20.
#' 
#' @export
#' 
#' @return A list, where each element is a PC, and each element is a named vector
#'         containing the top genes and their loadings.
#' 
extract_pc_loadings <- function(pca_obj,
                                pcs = 1:3,
                                top_n = 20) {
  # Check if pca_obj is a prcomp object
  if (!inherits(pca_obj, "prcomp")) {
    stop("pca_obj must be a prcomp object.")
  }

  # Check if pcs is a vector of positive integers
  if (!is.numeric(pcs) || any(pcs %% 1 != 0) || any(pcs <= 0)) {
    stop("pcs must be a vector of positive integers.")
  }
  
  # Check if top_n is a positive integer
  if (!is.numeric(top_n) || top_n %% 1 != 0 || top_n <= 0) {
    stop("top_n must be a positive integer.")
  }

  # Check if pcs are within the range of principal components
  # If not, keep only valid components
  allowed_pcs <- 1:ncol(pca_obj$x)
  invalid_pcs <- setdiff(pcs, allowed_pcs)
  if (length(invalid_pcs) > 0) {
    warning(paste("Principal components", invalid_pcs, "are invalid."))
    warning("Keeping only valid principal components.")
    pcs <- intersect(pcs, allowed_pcs)
  }

  # Extract loadings for each principal component
  loadings <- as.data.frame(pca_obj$rotation)

  if (top_n > nrow(loadings)) {
    warning("top_n is greater than the number of genes. Returning all genes.")
    top_n <- nrow(loadings)
  }

  # Get top genes for each principal component
  top_genes <- list()
  for (pc in pcs) {
    pc_loadings <- loadings[, pc]
    names(pc_loadings) <- rownames(loadings)

    # Sort loadings by absolute value
    top_genes[[paste0("PC_", pc)]] <- pc_loadings[order(abs(pc_loadings), decreasing = TRUE)][1:top_n]
  }

  return(top_genes)
}
