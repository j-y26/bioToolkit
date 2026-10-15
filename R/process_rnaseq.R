# Common functionalities for processing RNA-seq data

#' Columns of a DE table used by the DE table helpers
#'
#' Finds the log2 fold change and adjusted p-value columns of a DE table for
#' \code{rank_de_genes()}, \code{generate_de_rank_file()} and
#' \code{filter_de_genes()}, from the same presets as
#' \code{\link{standardize_de}}: columns given explicitly win; otherwise they
#' come from `method`, detected from the column names when "auto".
#'
#' @param data A DE table
#'
#' @param method "auto", "edgeR", "DESeq2", "limma" or "seurat" (case is
#'               ignored, as the old `de_method` argument did)
#'
#' @param logfc.col,padj.col Columns given by the user, or NULL
#'
#' @return A named character vector with `logfc` and `padj`
#'
#' @keywords internal
de_helper_columns <- function(data,
                              method,
                              logfc.col,
                              padj.col) {
    methods <- c("auto", names(DE_COLUMN_PRESETS))
    # The default (the whole vector) means "auto"; matching ignores case so
    # that the old de_method = "deseq2" keeps working
    method <- if (length(method) > 1) "auto" else methods[match(tolower(method), tolower(methods))]
    if (is.na(method)) {
        stop("`method` must be one of ", paste0("\"", methods, "\"", collapse = ", "), call. = FALSE)
    }
    if (!is.data.frame(data)) data <- as.data.frame(data)
    if (method == "auto" && !(!is.null(logfc.col) && !is.null(padj.col))) {
        # Only the two columns used here need to be present (a trimmed edgeR
        # table with logFC and FDR is enough, as before 0.10.0)
        found <- names(Filter(function(p) all(p[c("logfc", "padj")] %in% names(data)), DE_COLUMN_PRESETS))
        if (length(found) == 0) {
            stop("Could not detect the DE method from the columns of `data`; ",
                 "set `method` or `logfc.col` and `padj.col`", call. = FALSE)
        }
        method <- found[1]
        message("Detected ", method, " results (",
                paste(DE_COLUMN_PRESETS[[method]][c("logfc", "padj")], collapse = ", "), ")")
    }
    preset <- if (method == "auto") c(logfc = NA, padj = NA) else DE_COLUMN_PRESETS[[method]]
    cols <- c(logfc = logfc.col %||% unname(preset[["logfc"]]),
              padj = padj.col %||% unname(preset[["padj"]]))
    if (anyNA(cols)) stop("Set both `logfc.col` and `padj.col`", call. = FALSE)
    missing <- setdiff(cols, names(data))
    if (length(missing)) {
        stop("Column ", paste0("'", missing, "'", collapse = ", "), " not found in `data`", call. = FALSE)
    }
    cols
}


# Old (before 0.10.0) -> new argument names of the DE table helpers
DE_HELPER_ARG_RENAMES <- c(
    de_results = "data",
    de_method = "method",
    pvalue_col = "padj.col",
    log2fc_col = "logfc.col",
    log2fc_min = "logfc.cutoff",
    pvalue_max = "padj.cutoff",
    save_file = "file"
)


#' @title Calculate ranks of genes for differential expression analysis
#'
#' @description Ranks the genes of a DE table by a signed significance score,
#'              the pi-value:
#'              \deqn{pi = -log10(padj) * sign(logFC)}
#'              so that significantly up-regulated genes come first and
#'              significantly down-regulated genes last; a ranking for GSEA
#'              (see \code{\link{generate_de_rank_file}}).
#'
#' @param data A data frame of DE results with a log2 fold change and an
#'             adjusted p-value column (edgeR, DESeq2, limma or Seurat
#'             results, or any table with the columns given)
#'
#' @param method The tool that produced the results: "auto" (default) detects
#'               it from the column names, as \code{\link{standardize_de}}
#'               does; or one of "edgeR", "DESeq2", "limma" and "seurat"
#'
#' @param padj.col Column of the (adjusted) p-values used for the ranking. If
#'                 NULL (default), taken from `method` (e.g. "FDR" for edgeR,
#'                 "padj" for DESeq2)
#'
#' @param logfc.col Column of the log2 fold changes. If NULL (default), taken
#'                  from `method`
#'
#' @param ... Deprecated names of the arguments (`de_results`, `de_method`,
#'            `pvalue_col`, `log2fc_col`), accepted with a warning; renamed in
#'            0.10.0
#'
#' @details Ties are common among genes with the same adjusted p-value (e.g.
#'          many non-significant genes at padj = 1); they are broken by row
#'          order (`ties.method = "first"`).
#'
#' @export
#'
#' @return The input data frame, ordered by rank, with two added columns:
#'          \item{piValue}{The pi-value of the gene}
#'          \item{rank}{The rank of the gene based on the pi-value}
#'
#' @examples
#' de <- data.frame(
#'     logFC = c(2.1, -1.4, 0.2), PValue = c(1e-8, 1e-4, 0.6),
#'     FDR = c(3e-8, 2e-4, 0.6), row.names = c("ISG15", "CD3E", "ACTB")
#' )
#' rank_de_genes(de)
#'
rank_de_genes <- function(
    data,
    method = c("auto", "edgeR", "DESeq2", "limma", "seurat"),
    padj.col = NULL,
    logfc.col = NULL,
    ...) {
  # Old argument names (before 0.10.0) are mapped to the new ones
  list2env(renamed_args(list(...), DE_HELPER_ARG_RENAMES[1:4], "rank_de_genes",
                        names(match.call())), environment())
  cols <- de_helper_columns(data, method, logfc.col, padj.col)

  # Calculate pi-value
  data$piValue <- -log10(data[[cols[["padj"]]]]) * sign(data[[cols[["logfc"]]]])

  # Rank genes based on pi-value, from max to min
  data$rank <- rank(-data$piValue, ties.method = "first")

  # re-order the rows based on rank
  data <- data[order(data$rank), ]

  return(data)
}


#' @title Generate rank file for GSEA
#'
#' @description Ranks the genes of a DE table by their pi-value (see
#'              \code{\link{rank_de_genes}}) and returns, and optionally
#'              writes, the two-column table read by the GSEA desktop
#'              application (`.rnk`). For fgsea, use
#'              `setNames(rnk$rank, rnk$GeneName)`.
#'
#' @param file Path of the tab-separated rank file to write. Default is NULL
#'             (nothing written)
#'
#' @param ... Deprecated names of the arguments (`de_results`, `save_file`,
#'            `de_method`, `pvalue_col`, `log2fc_col`), accepted with a
#'            warning; renamed in 0.10.0
#'
#' @inheritParams rank_de_genes
#'
#' @export
#'
#' @return A data frame with columns `GeneName` (the row names of `data`) and
#'         `rank` (the pi-value), from highest to lowest; written to `file`
#'         when given
#'
#' @importFrom utils write.table
#' @import dplyr
#' @importFrom tibble rownames_to_column
#'
generate_de_rank_file <- function(
    data,
    file = NULL,
    method = c("auto", "edgeR", "DESeq2", "limma", "seurat"),
    padj.col = NULL,
    logfc.col = NULL,
    ...) {
  # Old argument names (before 0.10.0) are mapped to the new ones
  list2env(renamed_args(list(...), DE_HELPER_ARG_RENAMES[c(1:4, 7)], "generate_de_rank_file",
                        names(match.call())), environment())

  # Rank the genes based on differential expression results
  ranked_results <- rank_de_genes(data, method = method, padj.col = padj.col,
                                  logfc.col = logfc.col)

  # Select relevant columns for GSEA
  gsea_rank_file <- ranked_results %>%
    tibble::rownames_to_column(var = "GeneName") %>%
    mutate(rank = piValue) %>%
    select(GeneName, rank)

  # Save the rank file to disk if a file is given
  if (!is.null(file)) {
    write.table(gsea_rank_file, file = file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
  return(gsea_rank_file)
}


#' @title Split DE results into up- and down-regulated genes
#'
#' @description Keeps the genes of a DE table with an adjusted p-value below
#'              `padj.cutoff` and an absolute log2 fold change above
#'              `logfc.cutoff`, split by direction.
#'
#' @param logfc.cutoff Minimum absolute log2 fold change (strictly greater).
#'                     Default is 0 (any change)
#'
#' @param padj.cutoff Maximum adjusted p-value (strictly smaller). Default is
#'                    0.05
#'
#' @param ... Deprecated names of the arguments (`de_results`, `log2fc_min`,
#'            `pvalue_max`, `de_method`, `pvalue_col`, `log2fc_col`), accepted
#'            with a warning; renamed in 0.10.0
#'
#' @inheritParams rank_de_genes
#'
#' @export
#'
#' @return A list of two data frames, the rows of `data` that pass:
#'          \item{up_reg}{up-regulated genes}
#'          \item{down_reg}{down-regulated genes}
#'
#' @examples
#' de <- data.frame(
#'     logFC = c(2.1, -1.4, 0.2), PValue = c(1e-8, 1e-4, 0.6),
#'     FDR = c(3e-8, 2e-4, 0.6), row.names = c("ISG15", "CD3E", "ACTB")
#' )
#' filter_de_genes(de, logfc.cutoff = 1)
#'
filter_de_genes <- function(
    data,
    logfc.cutoff = 0,
    padj.cutoff = 0.05,
    method = c("auto", "edgeR", "DESeq2", "limma", "seurat"),
    padj.col = NULL,
    logfc.col = NULL,
    ...) {
  # Old argument names (before 0.10.0) are mapped to the new ones
  list2env(renamed_args(list(...), DE_HELPER_ARG_RENAMES[1:6], "filter_de_genes",
                        names(match.call())), environment())
  cols <- de_helper_columns(data, method, logfc.col, padj.col)
  logfc <- data[[cols[["logfc"]]]]
  padj <- data[[cols[["padj"]]]]

  # Filter genes based on log2 fold change and adjusted p-value; NA p-values
  # (e.g. DESeq2's independent filtering) never pass
  sig <- !is.na(padj) & padj < padj.cutoff
  up_reg <- data[sig & logfc > logfc.cutoff, ]
  down_reg <- data[sig & logfc < -logfc.cutoff, ]

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
#' @export
#'
#' @return A data frame containing the filtered raw counts
#'
filter_low_expression_genes <- function(raw_counts,
                                        min_group_size = 4) {
  check_installed("edgeR", reason = "for CPM-based gene filtering")

  # Convert raw counts to DGEList object
  dge <- edgeR::DGEList(counts = raw_counts)

  # Calculate counts per million (CPM)
  cpm_values <- edgeR::cpm(dge)

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
  allowed_pcs <- seq_len(ncol(pca_obj$x))
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
