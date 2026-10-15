# Gene set formats
#
# format_gene_sets() converts between the common ways gene sets are stored
# (GMT file, named list, GSEABase GeneSetCollection). It was kept from the
# former R/module_scoring_rnaseq.R when the module scoring functions were
# retired in 0.10.0 (backup in dev/legacy/).

#' @title Convert gene sets between formats
#'
#' @description Converts gene sets between a GMT file, a named list of gene
#'              vectors (the input of fgsea and of
#'              \code{Seurat::AddModuleScore()}) and a GSEABase
#'              `GeneSetCollection`.
#'
#' @param gsc A GeneSetCollection object or a list of gene sets or a path to a GMT file.
#'
#' @param to A character string indicating the format to convert to. Available
#'           options are "GeneSetCollection", "list", and "gmt".
#'
#' @param file For `to = "gmt"`, a file to write the GMT lines to. Default is
#'             NULL (only returned)
#'
#' @return A GeneSetCollection object, a named list of gene sets, or, for
#'         "gmt", the GMT lines (one per set: name, an empty description,
#'         then the genes, separated by tabs), invisibly when written to
#'         `file`
#'
#' @examples
#' sets <- list(SET_A = c("CD3E", "CD3D"), SET_B = c("MS4A1", "CD79A", "CD19"))
#' format_gene_sets(sets, to = "gmt")
#'
#' @export
#'
format_gene_sets <- function(gsc,
                             to = c("GeneSetCollection", "list", "gmt"),
                             file = NULL) {
  to <- match.arg(to)

  check_installed("GSEABase", reason = "for gene set manipulation")

  # Convert input to GeneSetCollection, if not already
  if (is.character(gsc)) {
    if (!file.exists(gsc)) {
      stop("The path to the GMT file does not exist")
    }
    gsc <- GSEABase::getGmt(gsc)
  } else if (is(gsc, "GeneSetCollection")) {
    # Do nothing

  } else if (is.list(gsc)) {
    gene_set_list <- lapply(names(gsc), function(set_name) {
      # A gene listed twice in a set (e.g. msigdbr symbols with several
      # Entrez IDs) is kept once: GeneSet() refuses duplicates
      genes <- unique(as.character(gsc[[set_name]]))
      GSEABase::GeneSet(genes, setName = set_name)
    })
    gsc <- GSEABase::GeneSetCollection(gene_set_list)
  } else {
    stop("gsc must be a GeneSetCollection object, a list of gene sets, or a path to a GMT file")
  }

  # Convert to desired format
  if (to == "GeneSetCollection") {
    return(gsc)
  } else if (to == "list") {
    gene_set_list <- lapply(names(gsc), function(set_name) {
      gene_set <- gsc[[set_name]]
      GSEABase::geneIds(gene_set)
    })
    names(gene_set_list) <- names(gsc)
    return(gene_set_list)
  } else if (to == "gmt") {
    # GMT lines built here: GSEABase::toGmt() without a file only prints them
    # (and returns NULL)
    gmt <- vapply(names(gsc), function(set_name) {
      paste(c(set_name, "", GSEABase::geneIds(gsc[[set_name]])), collapse = "\t")
    }, character(1), USE.NAMES = FALSE)
    if (!is.null(file)) {
      writeLines(gmt, file)
      return(invisible(gmt))
    }
    return(gmt)
  }
}
