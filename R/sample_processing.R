# Processing of the sample data

#' @title Extracting sample info based on regex keywords
#'
#' @description Extracting sample info based on regex keywords found in a
#'              string. If multiple keywords are found, the first keyword
#'              found will be used to extract the info.
#'
#' @param sample A string containing the sample info, must be a single string
#'
#' @param keywords A vector of regex keywords to search for
#'
#' @param info A vector of info represented by the keywords, must be the same
#'             length as keywords
#'
#' @return the info extracted from the sample based on the keywords found
#'
#' @export
#'
#' @examples
#' # Extracting sample info based on regex keywords
#' sample <- "Sample_2wks_1"
#' keywords <- c("[0-9]wks", "p7")
#' info <- c("postnatal_2wks", "postnatal_day_7")
#'
#' # Extracting the stage info
#' keywords_to_info(sample, keywords, info)
#'
keywords_to_info <- function(sample, keywords, info) {
  # Check if the length of keywords and info are the same
  if (length(keywords) != length(info)) {
    stop("The length of keywords and info must be the same")
  }

  # Check if any keywords are empty or NA or NULL
  if (any(is.na(keywords) | is.null(keywords) | keywords == "")) {
    stop("The keywords cannot be empty, NA or NULL")
  }

  # Check if the sample is a single string
  if (length(sample) != 1) {
    stop("The sample must be a single string")
  }

  # Extract the info based on the keywords
  for (i in seq_along(keywords)) {
    if (grepl(keywords[i], sample)) {
      return(info[i])
    }
  }
}

# Vectorized the function
keywords_to_info <- Vectorize(keywords_to_info, vectorize.args = "sample")
