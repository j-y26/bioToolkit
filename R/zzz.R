# Package infrastructure helpers

#' Check if a package is installed, with informative error message
#'
#' @param pkg Character string, the package name.
#' @param reason Optional character string explaining why the package is needed.
#'
#' @return Invisible TRUE if package is available; otherwise stops with an error.
#'
#' @keywords internal
check_installed <- function(pkg, reason = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- paste0("Package '", pkg, "' is required")
    if (!is.null(reason)) {
      msg <- paste0(msg, " ", reason)
    }
    # Determine install command
    bioc_pkgs <- c(
      "edgeR", "Matrix", "GSEABase", "ChIPseeker", "GenomicRanges",
      "GSVA", "BiocParallel", "IRanges", "Rsamtools", "rtracklayer",
      "Biobase"
    )
    if (pkg %in% bioc_pkgs) {
      msg <- paste0(
        msg,
        ".\nInstall it with: BiocManager::install('", pkg, "')"
      )
    } else {
      msg <- paste0(msg, ".\nInstall it with: install.packages('", pkg, "')")
    }
    stop(msg, call. = FALSE)
  }
  invisible(TRUE)
}


#' Check if an object is a ggplot
#'
#' @param x An object to check.
#'
#' @return Logical, TRUE if the object inherits from "gg".
#'
#' @keywords internal
is_ggplot <- function(x) {
  inherits(x, "gg")
}
