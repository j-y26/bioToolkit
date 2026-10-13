# Source data behind figures (roadmap item F6)
#
# Design
# ------
# Journals ask for the numbers behind each figure panel. Every plot function
# of the package keeps the tidy table it built: either as the plot's own data
# (`plot$data`), or, when its layers use different data, as the
# "source_data" attribute, which takes precedence. `export_source_data()`
# writes that table out, one table per plot of a patchwork.


#' @title Export the data behind a plot
#'
#' @description Writes the table each plot of the package is built from, so
#'              that it can be submitted as source data with a figure. For a
#'              patchwork, each plot is written separately: one CSV file per
#'              plot, or one sheet per plot of an XLSX file.
#'
#'              The table is the "source_data" attribute of a plot when it has
#'              one, else its data (`plot$data`); row names (e.g. cell names)
#'              are kept as a first column.
#'
#' @param plot A ggplot or patchwork object
#'
#' @param file Output file, ".csv" or ".xlsx". For a patchwork written to CSV,
#'             the plot names (or numbers) are appended to the file name. If
#'             NULL (default), nothing is written and the tables are only
#'             returned
#'
#' @return A named list of data frames, one per plot, invisibly
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
#' f <- tempfile(fileext = ".csv")
#' export_source_data(p, f)
#' head(read.csv(f))
#'
#' @export
#'
export_source_data <- function(plot, file = NULL) {
    tables <- source_tables(plot)
    if (length(tables) == 0) stop("The plot holds no data to export", call. = FALSE)
    if (is.null(file)) return(invisible(tables))

    ext <- tolower(sub(".*\\.", "", basename(file)))
    if (ext == "csv") {
        files <- if (length(tables) == 1) {
            file
        } else {
            paste0(sub("\\.[^.]*$", "", file), "_", safe_file_part(names(tables)), ".csv")
        }
        for (i in seq_along(tables)) utils::write.csv(tables[[i]], files[i], row.names = FALSE)
    } else if (ext == "xlsx") {
        check_installed("writexl", reason = "to write XLSX files")
        # Sheet names are limited to 31 characters and must be unique
        sheets <- make.unique(substr(names(tables), 1, 28), sep = "_")
        writexl::write_xlsx(stats::setNames(tables, sheets), file)
    } else {
        stop("`file` must end in .csv or .xlsx", call. = FALSE)
    }
    invisible(tables)
}


#' The source tables of a plot, one per ggplot of a patchwork
#'
#' @param plot A ggplot or patchwork object
#'
#' @param prefix Name of the plot, for nested patchworks
#'
#' @return A named list of data frames
#'
#' @keywords internal
source_tables <- function(plot, prefix = NULL) {
    if (inherits(plot, "patchwork")) {
        plots <- patchwork_parts(plot)$plots
        nm <- names(plots)
        if (is.null(nm)) nm <- rep("", length(plots))
        nm[!nzchar(nm)] <- paste0("panel", seq_along(plots))[!nzchar(nm)]
        if (!is.null(prefix)) nm <- paste(prefix, nm, sep = "_")
        out <- unlist(lapply(seq_along(plots), function(i) source_tables(plots[[i]], nm[i])),
                      recursive = FALSE)
        return(out)
    }
    if (!is_ggplot(plot)) return(list())
    data <- attr(plot, "source_data")
    if (is.null(data)) data <- plot$data
    # A plot with no data of its own (e.g. layers each with their own data, or
    # a spacer) contributes nothing
    if (!is.data.frame(data) || nrow(data) == 0) return(list())
    data <- as.data.frame(data)
    if (.row_names_info(data) > 0) {
        data <- cbind(data.frame(row = rownames(data), stringsAsFactors = FALSE), data)
        rownames(data) <- NULL
    }
    stats::setNames(list(data), if (is.null(prefix)) "plot" else prefix)
}


#' Make strings safe as parts of file names
#'
#' @param x Strings
#'
#' @keywords internal
safe_file_part <- function(x) {
    gsub("[^A-Za-z0-9._-]+", "_", x)
}
