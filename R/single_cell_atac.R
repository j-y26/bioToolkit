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
#' @param palette A valid palette name from "alphabet", "alphabet2", "glasbey", "polychrome", "stepped", and "parade".
#'                Alternatively, a vector of colors can be provided with the same length as unique groups.
#'
#' @param plot_title A character string for the title of the plot. Default is
#'                   "Peak Category QC".
#'
#' @return A ggplot object representing the bar plot of peak categories.
#'
#' @import ggplot2 patchwork
#' @importFrom dplyr bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom ChIPseeker annotatePeak
#'
#' @export
#'
peak_category_plot <- function(peak_list,
                               peak_names = NULL,
                               TxDb,
                               palette = "polychrome",
                               plot_title = "Peak Category QC") {
  # Annotate peaks
  annotated_peaks <- peak_category_qc(peak_list, peak_names, TxDb)

  # Combine annotations into a single data frame
  df_list <- lapply(annotated_peaks, function(x) {
    df <- as.data.frame(x@annoStat)
    df$count <- df$Frequency / 100 * x@peakNum
    return(df)
  })

  df_combined <- dplyr::bind_rows(df_list, .id = "peak_set")

  # Plot the data
  p_stack <- grouped_bar_plot(
    df = df_combined,
    group.by = "Feature",
    split.by = "peak_set",
    frequency = "count",
    mode = "stack",
    palette = palette,
    show_totals = TRUE,
    title = plot_title
  ) +
    theme(legend.position = "None")

  p_fill <- grouped_bar_plot(
    df = df_combined,
    group.by = "Feature",
    split.by = "peak_set",
    frequency = "count",
    mode = "fill",
    palette = palette,
    show_totals = FALSE
  )

  return(p_stack + p_fill)
}


#' @title Peak TSS enrichment
#'
#' @description This function calculates the TSS (Transcription Start Site) enrichment
#' for a given set of peaks. It uses the `ChIPseeker` package to annotate the peaks
#' and then computes the TSS enrichment.
#'
#' @param peaks A list of GRanges objects, each representing peaks identified
#'              by a different peak calling method or sample.
#'
#' @param peak_names A character vector of names corresponding to the peaks in
#'                  `peaks`. If not provided, names will be extracted from
#' '                `peaks` names. Default is `NULL`.
#'
#' @param TxDb A TxDb object containing gene annotations.
#'
#' @param upstream A numeric value specifying the upstream region to consider for TSS enrichment.
#'
#' @param downstream A numeric value specifying the downstream region to consider for TSS enrichment.
#'
#' @param plot_title A character string for the title of the plot. Default is
#'                   "Peak TSS Enrichment".
#'
#' @return A ggplot object representing the TSS enrichment plot.
#'
#' @import ggplot2 patchwork
#' @importFrom aplot plot_list
#' @importFrom ChIPseeker peak_Profile_Heatmap
#'
#' @export
#'
peak_tss_enrichment <- function(peaks,
                                peak_names = NULL,
                                TxDb,
                                upstream = 3000,
                                downstream = 3000,
                                plot_title = "Peak TSS Enrichment") {
  # Check if peaks is a named list
  if (!is.list(peaks)) {
    stop("peaks must be a list of GRanges objects")
  }
  if (is.null(names(peaks)) && is.null(peak_names)) {
    stop("peaks must be a named list or peak_names must be provided")
  }
  if (!is.null(peak_names) && length(peak_names) != length(peaks)) {
    stop("peak_names must be the same length as peaks")
  }
  if (!is.null(peak_names)) {
    names(peaks) <- peak_names
  }

  # Check if GRanges objects are in the list
  if (!all(sapply(peaks, inherits, "GRanges"))) {
    stop("All elements in peaks must be GRanges objects")
  }

  # Check if TxDb is provided and is a valid TxDb object
  if (missing(TxDb) || !inherits(TxDb, "TxDb")) {
    stop("TxDb must be provided and must be a valid TxDb object")
  }

  # Calculate TSS enrichment
  p_list <- lapply(seq_along(peaks), function(i) {
    peak <- peaks[[i]]
    name <- names(peaks)[i]
    cat("Processing peak set:", name, "\n")

    p <- ChIPseeker::peak_Profile_Heatmap(
      peak,
      TxDb = TxDb,
      upstream = upstream,
      downstream = downstream,
      title = name,
      xlab = "Distance to TSS (bp)",
    )
    return(p)
  })

  # Combine plots into a single plot
  p_combined <- aplot::plot_list(
    gglist = p_list,
    nrow = 1
  ) +
    patchwork::plot_annotation(plot_title)
  return(p_combined)
}


#' @title Peak width distribution
#'
#' @description This function calculates the distribution of peak widths for a
#'              given set of peaks. It generates a histogram with a density
#'              curve overlay to visualize the distribution of peak widths. A
#'              boxplot is also generated to show the distribution of peak
#'              widths across different peak sets.
#'
#' @param peaks A list of GRanges objects, each representing peaks identified
#'              by a different peak calling method or sample.
#'
#' @param peak_names A character vector of names corresponding to the peaks in
#'                   `peaks`. If not provided, names will be extracted from
#'                   `peaks` names. Default is `NULL`.
#'
#' @param plot_title A character string for the title of the plot. Default is
#                   "Peak Width Distribution".
#'
#' @return A patchwork object containing a histogram of peak widths with a density
#'         curve overlay and a boxplot showing the distribution of peak widths
#'         across different peak sets.
#'
#' @import ggplot2 patchwork
#' @importFrom GenomicRanges width
#'
#' @export
#'
peak_width_distribution <- function(peaks,
                                    peak_names = NULL,
                                    plot_title = "Peak Width Distribution") {
  # Check if peaks is a named list
  if (!is.list(peaks)) {
    stop("peaks must be a list of GRanges objects")
  }
  if (is.null(names(peaks)) && is.null(peak_names)) {
    stop("peaks must be a named list or peak_names must be provided")
  }
  if (!is.null(peak_names) && length(peak_names) != length(peaks)) {
    stop("peak_names must be the same length as peaks")
  }
  if (!is.null(peak_names)) {
    names(peaks) <- peak_names
  }

  # Check if GRanges objects are in the list
  if (!all(sapply(peaks, inherits, "GRanges"))) {
    stop("All elements in peaks must be GRanges objects")
  }

  cols <- ggplotColours(length(peaks))

  # Determine the max and min widths to scale the plots
  max_width <- max(sapply(peaks, function(x) max(GenomicRanges::width(x), na.rm = TRUE)))
  min_width <- min(sapply(peaks, function(x) min(GenomicRanges::width(x), na.rm = TRUE)))
  # Round to the nearest 50
  max_width <- ceiling(max_width / 50) * 50
  min_width <- floor(min_width / 50) * 50

  # Plot for each peak set
  p_list <- lapply(seq_along(peaks), function(i) {
    peak <- peaks[[i]]
    name <- names(peaks)[i]

    # Create a data frame for plotting
    df <- data.frame(
      width = GenomicRanges::width(peak),
      peak_set = name
    )

    # Histogram with density curve
    p <- ggplot(df, aes(x = width)) +
      geom_histogram(aes(y = ..density..), binwidth = 50, fill = cols[i], color = "black", alpha = 0.6, boundary = 0) +
      geom_density(color = "darkblue", size = 1) +
      labs(title = name,
           x = "Peak Width (bp)",
           y = "Density") +
      theme_classic()

    # Overlay boxplot
    plot_data <- ggplot_build(p)
    y_hist <- plot_data$data[[1]]$y
    y_dens <- plot_data$data[[2]]$y
    y_max <- max(c(y_hist, y_dens), na.rm = TRUE)

    p <- p +
      geom_boxplot(aes(y = -0.08 * y_max), width = 0.1 * y_max, fill = cols[i], outlier.shape = NA, alpha = 0.6) +
      scale_x_continuous(limits = c(min_width, max_width))

    return(p)
  })

  # Combine plots into a single plot
  p_combined <- patchwork::wrap_plots(p_list, ncol = 1, axes = "collect") +
    patchwork::plot_annotation(title = plot_title)

  return(p_combined)
}


#' @title Calculate FRiP score for peaks
#'
#' @description This function calculates the FRiP (Fraction of Reads in Peaks) score
#' for a given set of peaks. It computes the fraction of reads that overlap with
#' the peaks and returns a data frame with the FRiP scores for each peak set.
#'
#' @param peaks A list of GRanges objects, each representing peaks identified
#'              by a different peak calling method or sample.
#'
#' @param peak_names A character vector of names corresponding to the peaks in
#'                   `peaks`. If not provided, names will be extracted from
#'                  `peaks` names. Default is `NULL`.
#'
#' @param fragment_file A character string specifying the path to the fragment file
#'                      (in BED format) containing the reads. This file should
#'                      be deduplicated and sorted. This file can be either
#'                      Tabix indexed or not. However, for large files, it is
#'                      recommended to use a Tabix indexed file for faster access.
#'
#' @param threads An integer specifying the number of threads to use for parallel processing.
#'                Default is 8.
#'
#' @param is_indexed A logical value indicating whether the fragment file is Tabix indexed.
#'                   Default is to detect automatically (NULL).
#'
#' @return A data frame containing the FRiP scores for each peak set.
#'
#' @importFrom GenomicRanges GRanges seqnames sort
#' @importFrom IRanges IRanges
#' @importFrom Rsamtools TabixFile scanTabix
#' @importFrom rtracklayer import
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom data.table fread
#'
#' @export
#'
peak_calc_frip <- function(peaks,
                           peak_names = NULL,
                           fragment_file,
                           threads = 8,
                           is_indexed = NULL) {
  # Check if peaks is a named list
  if (!is.list(peaks)) {
    stop("peaks must be a list of GRanges objects")
  }
  if (is.null(names(peaks)) && is.null(peak_names)) {
    stop("peaks must be a named list or peak_names must be provided")
  }
  if (!is.null(peak_names) && length(peak_names) != length(peaks)) {
    stop("peak_names must be the same length as peaks")
  }
  if (!is.null(peak_names)) {
    names(peaks) <- peak_names
  }

  if (!file.exists(fragment_file)) stop("Fragment file does not exist")

  # Check if GRanges objects are in the list
  if (!all(sapply(peaks, inherits, "GRanges"))) {
    stop("All elements in peaks must be GRanges objects")
  }

  # Detect index if not provided
  if (is.null(is_indexed)) {
    is_indexed <- file.exists(paste0(fragment_file, ".tbi"))
  }

  # Ensure peaks are sorted and have seqinfo
  peaks <- lapply(peaks, function(x) {
    x <- GenomicRanges::sort(x)
  })
  seqs <- unique(unlist(lapply(peaks, function(p) as.character(GenomicRanges::seqnames(p)))))

  # Set up parallel processing
  bp_param <- BiocParallel::MulticoreParam(workers = threads, progressbar = TRUE)

  # Pre-index peaks per chromosome
  peaks_by_chr <- lapply(seqs, function(chr) {
    chr_peaks <- lapply(peaks, function(x) x[GenomicRanges::seqnames(x) == chr])
    names(chr_peaks) <- names(peaks)
    return(chr_peaks)
  })
  names(peaks_by_chr) <- seqs

  # Define a function to calculate FRiP for a single chromosome
  process_chr <- function(chr) {
    chr_peaks <- peaks_by_chr[[chr]]
    cat("Processing chromosome:", chr, "in process", Sys.getpid(), "\n")

    # Read the fragment file for the current chromosome depending on whether it is indexed or not
    if (is_indexed) {
      tabix_file <- Rsamtools::TabixFile(fragment_file)
      param <- GenomicRanges::GRanges(chr, IRanges::IRanges(1, MAX_TABIX_END)) # Max Tabix end is beyond standard chromosome lengths
      cat("Tabix file detected. Using Tabix to read fragments for chromosome")
      cat("Reading fragments:", chr, "\n")
      fragments <- Rsamtools::scanTabix(tabix_file, param = param)[[1]]
      if (length(fragments) == 0) {
        warning(paste("No fragments found for chromosome", chr, "in the fragment file."))
          empty_result <- rep(list(c(total = 0L, in_peaks = 0L)), length(chr_peaks))
          setNames(empty_result) <- names(chr_peaks)
        return(empty_result)
      }
      dt <- data.table::fread(text = fragments, header = FALSE, sep = "\t", showProgress = TRUE)
      # Convert to GRanges
      dt <- GenomicRanges::GRanges(
        seqnames = dt$V1,
        ranges = IRanges::IRanges(start = dt$V2 + 1, end = dt$V3),
        name = dt$V4,
        score = dt$V5
      )
      cat("Fragments read successfully for chromosome:", chr, "\n")
    } else {
      cat("No Tabix index detected. Reading fragments directly from the file.\n")
      cat("Reading fragments:", chr, "\n")
      dt <- rtracklayer::import(
        fragment_file,
        which = GenomicRanges::GRanges(chr, IRanges::IRanges(1, MAX_TABIX_END)),
        format = "BED"
      )
      if (is.null(td) || length(td) == 0) {
        warning(paste("No fragments found for chromosome", chr, "in the fragment file."))
        empty_result <- rep(list(c(total = 0L, in_peaks = 0L)), length(chr_peaks))
        setNames(empty_result) <- names(chr_peaks)
        return(empty_result)
      }
      cat("Fragments read successfully for chromosome:", chr, "\n")
    }

    # Calculate overlaps with peaks
    cat("Calculating overlaps:", chr, "\n")
    total_frag_chr <- length(dt)
    in_peaks_chr <- sapply(chr_peaks, function(peak) {
      if (length(peak) == 0) return(0L)
      overlaps <- GenomicRanges::countOverlaps(dt, peak, type = "any")
      return(sum(overlaps > 0))
    })

    return(mapply(function(in_p) c(total = total_frag_chr, in_peaks = in_p), in_peaks_chr, SIMPLIFY = FALSE))
  }

  # Process each chromosome in parallel
  results <- BiocParallel::bplapply(seqs, process_chr, BPPARAM = bp_param)

  warning(paste0("An maximum length of ", MAX_TABIX_END, " is set for reading the
                   fragment file. This value is beyond standard chromosome
                   length of eukaryotic orgnisms. You can safely ignore this
                   warning if this value is beyond the length of the chromosomes
                   in your organism."))

  # Aggregate results
  total_all <- rep(0L, length(peaks))
  in_peaks_all <- rep(0L, length(peaks))

  for (chr_res in results) {
    for (i in seq_along(chr_res)) {
      total_all[i] <- total_all[i] + chr_res[[i]]["total"]
      in_peaks_all[i] <- in_peaks_all[i] + chr_res[[i]]["in_peaks"]
    }
  }

  # Combine results into a data frame
  frip_df <- data.frame(
    peak_set = names(peaks),
    total_fragments = total_all,
    fragments_in_peaks = in_peaks_all,
    frip_score = ifelse(total_all > 0, in_peaks_all / total_all, NA_real_)
  )

  return(frip_df)
}


#' @title Plot FRiP scores
#' 
#' @description This function generates a bar plot of FRiP (Fraction of Reads in Peaks)
#' 
#' @param frip_df A data frame containing the FRiP scores for each peak set, the output of `peak_calc_frip`.
#' 
#' @param palette A valid palette name from "alphabet", "alphabet2", "glasbey", "polychrome", "stepped", and "parade".
#' Alternatively, a vector of colors can be provided with the same length as unique groups.
#' 
#' @param plot_title A character string for the title of the plot. Default is "FRiP Score".
#' 
#' @return A ggplot object representing the bar plot of FRiP scores.
#' 
#' @import ggplot2
#' @importFrom Seurat DiscretePalette
#' 
#' @export
#' 
peak_plot_frip <- function(frip_df,
                           palette = "polychrome",
                           plot_title = "FRiP Score") {
  # Check if frip_df is a data frame
  if (!is.data.frame(frip_df)) {
    stop("frip_df must be a data frame")
  }

  # Check required columns
  required_cols <- c("peak_set", "frip_score")
  if (!all(required_cols %in% colnames(frip_df))) {
    stop(paste("frip_df must contain the following columns:", paste(required_cols, collapse = ", ")))
  }

  # Create the bar plot
  p <- ggplot(frip_df, aes(x = peak_set, y = frip_score, fill = peak_set)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = Seurat::DiscretePalette(
      palette = palette,
      n = length(unique(frip_df$peak_set))
    )) +
    labs(title = plot_title, x = "Peak Set", y = "FRiP Score") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}