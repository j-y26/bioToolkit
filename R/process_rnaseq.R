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
