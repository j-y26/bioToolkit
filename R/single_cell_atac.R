# Functions for single-cell ATAC-seq analysis

#' @title Peak category QC across one or more peak calling methods/samples
#'
#' @description This function performs quality control on peak categories across
#'              one or more peak calling methods. It generates a bar plot
#'              showing the distribution of peaks across different categories.
#'
#' @param peak_list A named list of GRanges objects, each representing peaks
#'                  identified by a different peak calling method or sample.
#'
#' @param peak_names A character vector of names corresponding to the peaks in
#'                   `peak_list`. If not provided, names will be extracted from
#'                   `peak_list` names. Default is `NULL`.
#'
#' @param TxDb A TxDb object containing gene annotations. This is used to
#'             categorize peaks into promoter, intronic, intergenic, and
#'             exonic regions.
#'
#' @return A list of csAnno objects, each containing the detailed annotation of
#'         each peak set. See `ChIPseeker` package for more details.
#'
#' @import GenomicRanges
#' @importFrom ChIPseeker annotatePeak
#'
#' @export
#'
peak_category_qc <- function(peak_list,
                             peak_names = NULL,
                             TxDb) {
  # Check if peak_list is a named list
  if (!is.list(peak_list)) {
    stop("peak_list must be a list of GRanges objects")
  }
  if (is.null(names(peak_list)) && is.null(peak_names)) {
    stop("peak_list must be a named list or peak_names must be provided")
  }
  if (!is.null(peak_names) && length(peak_names) != length(peak_list)) {
    stop("peak_names must be the same length as peak_list")
  }
  if (!is.null(peak_names)) {
    names(peak_list) <- peak_names
  }

  # Check if GRanges objects are in the list
  if (!all(sapply(peak_list, inherits, "GRanges"))) {
    stop("All elements in peak_list must be GRanges objects")
  }

  # Check if TxDb is provided and is a valid TxDb object
  if (missing(TxDb) || !inherits(TxDb, "TxDb")) {
    stop("TxDb must be provided and must be a valid TxDb object")
  }

  # Annotate peaks
  annotated_peaks <- lapply(peak_list, function(peaks) {
    ChIPseeker::annotatePeak(peaks, TxDb = TxDb, tssRegion = c(-3000, 3000),
                 verbose = FALSE)
  })
  names(annotated_peaks) <- names(peak_list)

  return(annotated_peaks)
}


#' @title Plotting peak category across one or more peak calling methods/samples
#'
#' @description This function performs quality control on peak categories across
#'              one or more peak calling methods. It generates a bar plot
#'              showing the distribution of peaks across different categories.
#'
#' @param peak_list A named list of GRanges objects, each representing peaks
#'                  identified by a different peak calling method or sample.
#'
#' @param peak_names A character vector of names corresponding to the peaks in
#'                   `peak_list`. If not provided, names will be extracted from
#'                   `peak_list` names. Default is `NULL`.
#'
#' @param TxDb A TxDb object containing gene annotations. This is used to
#'             categorize peaks into promoter, intronic, intergenic, and
#'             exonic regions.
#'
#' @param plot_title A character string for the title of the plot. Default is
#'                   "Peak Category QC".
#'
#' @palette A character string specifying the color palette to use for the plot.
#'                Default is "Set1" from RColorBrewer.
#'
#' @return A ggplot object representing the bar plot of peak categories.
#'
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom ChIPseeker annotatePeak
#'
#' @export
#'
peak_category_plot <- function(peak_list,
                               peak_names = NULL,
                               TxDb,
                               plot_title = "Peak Category QC",
                               palette = "Set1") {
  # Annotate peaks
  annotated_peaks <- peak_category_qc(peak_list, peak_names, TxDb)

  # Combine annotations into a single data frame
  df_list <- lapply(annotated_peaks, function(x) {
    as.data.frame(x@annoStat)
  })
  combined_df <- dplyr::bind_rows(df_list, .id = "peak_set")

  # Create the bar plot
  p <- ggplot(combined_df, aes(x = peak_set, y = Frequency, fill = category)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = plot_title, x = "Peak Set", y = "Count") +
    scale_fill_brewer(palette = palette) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}
