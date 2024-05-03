# Common functionalities for processing RNA-seq data

#' @title Calculate ranks of genes for differential expression analysis
#' 
#' @description This function calculates ranks of genes based on the
#'              differential expression analysis results for RNA-seq data.
#' 
#' @param de_results A data frame containing the results of differential
#'                   expression analysis. The data frame should have at least
#'                   a column for log2 fold change and a column for (adjusted)
#'                   p-value.
#' 
#' @param pvalue_col The name of the column in \code{de_results} that contains
#'                   the p-values.
#' 
#' @param log2fc_col The name of the column in \code{de_results} that contains
#'                   the log2 fold changes.
#' 
#' 