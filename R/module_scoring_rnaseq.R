# Perform Gene Module Scoring on RNA-seq data

#' @title Calculate gene module score
#'
#' @description Calculate gene module score for each sample using
#'              average expression of genes in the module with
#'              background correction. This function is inspired by
#'              Seurat's \code{AddModuleScore} function and is adapted
#'              to be used with RNA-seq data.
#'
#' @param data A matrix of gene expression data with genes in rows
#'             and samples in columns. Rownames should be genes and column
#'             names should be samples. Must be a numeric matrix.
#'
#' @param modules A list of gene modules. Each element of the list should
#'                be a character vector of gene names. The names of the list
#'                should be the module names.
#'
#' @param normalized A logical value indicating whether the data is
#'                   normalized. If \code{TRUE}, the function will assume
#'                   that the data is already normalized by CPM or TPM or FPKM.
#'
#' @param pool A character vector of gene names to be used as background
#'             for the gene module scoring. If \code{NULL}, all genes in
#'             the data will be used as background.
#'
#' @param nbin Number of bins to divide genes into based on expression levels.
#'             Default is 24.
#'
#' @param ctrl Number of control genes to sample per gene in the gene set.
#'
#' @param seed Seed for reproducibility. Default is 123.
#'
#' @return A data frame with rows as samples and columns as module scores.
#'
#' @export
#'
#' @importFrom ggplot2 cut_number
#' @importFrom stats rnorm
#' @importFrom Matrix rowMeans colMeans
#' @importFrom Seurat CaseMatch
#' @importFrom edgeR DGEList calcNormFactors cpm
#'
calc_module_scores <- function(data,
                               modules,
                               normalized = TRUE,
                               pool = NULL,
                               nbin = 24,
                               ctrl = 100,
                               seed = 123) {

  # Check if data is a matrix
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("data must be a matrix or dataframe")
  }

  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }

  # Check if data is numeric
  if (!is.numeric(data)) {
    stop("data must be a numeric matrix")
  }

  # Check if normalized is a logical
  if (!is.logical(normalized)) {
    stop("normalized must be a logical value")
  }

  # Check if modules is a list
  if (!is.list(modules)) {
    stop("modules must be a list")
  }

  # Check if pool is a character vector
  if (!is.null(pool) && !is.character(pool)) {
    stop("pool must be a character vector")
  }

  # Check if nbin is numeric
  if (!is.numeric(nbin)) {
    stop("nbin must be a numeric value")
  }

  # Check if ctrl is a numeric
  if (!is.numeric(ctrl)) {
    stop("ctrl must be a numeric value")
  }

  # Check if seed is a numeric
  if (!is.numeric(seed)) {
    stop("seed must be a numeric value")
  }

  # Set seed
  set.seed(seed)

  # Get gene names
  genes <- rownames(data)

  # Get sample names
  samples <- colnames(data)

  # Get background genes
  pool <- pool %||% genes

  # Get gene expression
  if (!normalized) {
    # Assume using raw counts
    dge <- edgeR::DGEList(counts = data)
    dge <- edgeR::calcNormFactors(dge)
    data <- edgeR::cpm(dge, log = TRUE)
    cat("Raw data is normalized and log2-transformed\n")
  } else {
    # Assume using normalized data, log-transform it
    data <- log2(data + 1)
    cat("Normalized data is log2-transformed\n")
  }

  # Validate each module
  for (i in seq_along(modules)) {
    module <- modules[[i]]
    module <- Seurat::CaseMatch(module, genes)
    cat(paste0("Validating module ", names(modules)[i], "...\n"))
    module_diff <- setdiff(module, genes)
    if (length(module_diff) > 0) {
      warning(paste("- The following genes are not present in the data and will be ignored:", paste(module_diff, collapse = ", ")))
    }

    # Remove genes not in the data
    modules[[i]] <- intersect(module, genes)
  }
  cat("All modules are validated\n")

  # Checking that each module has at least 2 genes
  checked_length <- sapply(modules, function(x) length(x) >= 2)
  if (!all(checked_length)) {
    stop(paste("The following modules have less than 2 genes:",
               names(modules)[!checked_length],
                "Please remove these modules and try again.\n",
                "Exiting..."))
  }

  # Bin genes by expression levels
  pool_avg <- Matrix::rowMeans(data[pool, , drop = FALSE])
  pool_avg <- pool_avg[order(pool_avg)]
  pool_cut <- ggplot2::cut_number(pool_avg + stats::rnorm(length(pool_avg))/1e30,
                                  n = nbin,
                                  labels = FALSE,
                                  right = FALSE)
  names(pool_cut) <- names(pool_avg)

  # Define a list of control genes to use for each gene in the gene set
  control_genes <- vector(mode = "list", length = length(modules))
  for (i in seq_along(modules)) {
    module <- modules[[i]]
    for (j in seq_along(module)) {
      gene <- module[j]

      # TryCatch control gene sampling
      tryCatch({
        control_genes[[i]] <- c(control_genes[[i]],
                              names(sample(x = pool_cut[which(pool_cut == pool_cut[module[j]])],
                                           size = ctrl,
                                           replace = FALSE)))
      }, error = function(e) {
        warning(e$message)
        stop(paste("Error in sampling control genes for gene", gene, "in module", names(modules)[i], "\n",
                   "Reduce nbin or ctrl and try again. Exiting..."))
      })
    }
  }

  control_genes <- lapply(control_genes, unique)

  # Calculate control scores
  control_scores <- matrix(
    data = numeric(1L),
    nrow = length(control_genes),
    ncol = ncol(data)
  )

  for (i in seq_along(control_genes)) {
    control_scores[i, ] <- Matrix::colMeans(data[control_genes[[i]], , drop = FALSE])
  }

  # Calculate module scores
  module_scores <- matrix(
    data = numeric(1L),
    nrow = length(modules),
    ncol = ncol(data)
  )

  for (i in seq_along(modules)) {
    module_scores[i, ] <- Matrix::colMeans(data[modules[[i]], , drop = FALSE])
  }

  # Background correction
  module_scores_corr <- module_scores - control_scores

  # Create a data frame to store module scores
  module_scores_df <- data.frame(t(module_scores_corr))
  rownames(module_scores_df) <- samples
  colnames(module_scores_df) <- names(modules)

  set.seed(NULL)

  return(module_scores_df)
}


#' @title Calculate Gene module enrichment using ssGSEA
