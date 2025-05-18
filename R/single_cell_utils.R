# Utilities that are useful for single cell analysis

#' @title Generate a pivot table from cell metadata
#'
#' @description This function generates a pivot table from cell metadata and
#'              returns a data frame with the two cell metadata as rows and
#'              columns, with the values being either the count of cells or
#'              the percentage of cells (by rows or columns or total).
#'
#' @param row_metadata A vector of cell metadata to be used as rows, often a
#'                     metadata that represents cell clusters or cell types
#'
#' @param col_metadata A vector of cell metadata to be used as columns, must be
#'                     the same length as row_metadata
#'
#' @param type The type of pivot table to generate, either "percent" or "count"
#'
#' @param margin The margin to calculate the percentage, either "row", "column"
#'               or "total", used only when type = "percent"
#'
#' @return a data frame with the pivot table
#'
#' @import dplyr
#' @importFrom tidyr pivot_wider
#'
#' @export
#'
get_meta_pivot <- function(
  row_metadata,
  col_metadata,
  type = c("percent", "count"),
  margin = c("row", "column", "total")
) {
  # Check if the row_metadata and col_metadata are the same length
  if (length(row_metadata) != length(col_metadata)) {
    stop("The row_metadata and col_metadata must be the same length")
  }

  # Obtain settings
  type <- match.arg(type)
  margin <- match.arg(margin)

  # Create a count table
  count_table <- as.data.frame(table(row_metadata, col_metadata))

  # If the type is count, return the count table
  if (type == "count") {
    # Generate the count table
    count_table <- tidyr::pivot_wider(
      data = count_table,
      names_from = col_metadata,
      values_from = Freq
    )
    return(count_table)
  } else {
    # Generate the percentage table
    if (margin == "row") {
      pct_table <- count_table %>%
        group_by(row_metadata) %>%
        mutate(percent = Freq / sum(Freq) * 100) %>%
        ungroup()
    } else if (margin == "column") {
      pct_table <- count_table %>%
        group_by(col_metadata) %>%
        mutate(percent = Freq / sum(Freq) * 100) %>%
        ungroup()
    } else {
      pct_table <- count_table %>%
        mutate(percent = Freq / sum(Freq) * 100)
    }

    # Generate the percentage table
    pct_table <- pct_table %>%
      select(-Freq) %>%
      tidyr::pivot_wider(
        names_from = col_metadata,
        values_from = percent
      )
    return(pct_table)
  }
}


#' @title Tranforming a dimensionality reduction embedding, flipping and/or rotating
#'
#' @description This function transforms a dimensionality reduction embedding by
#'              flipping and/or rotating the embedding. This is useful for
#'              visualizing the embedding in a consistent way, for example, when
#'              comparing two embeddings or when the embedding is not in the
#'              desired orientation. This function will only work for 2D
#'              embeddings. It will  modify/overwrite the original embedding
#'              matrix. However, the biological meaning of the embedding is
#'              not changed. Only for visualization purposes.
#' 
#' @param obj A Seurat object
#' 
#' @param reduction The name of the dimensionality reduction to be transformed
#'                 (e.g. "pca", "umap", "tsne")
#' 
#' @param flip_horizontal A boolean indicating whether to flip the embedding
#'                        horizontally (default = FALSE)
#' 
#' @param flip_vertical A boolean indicating whether to flip the embedding
#'                       vertically (default = FALSE)
#' 
#' @param rotate A numeric value indicating the angle (in degrees) to rotate
#'               the embedding (default = 0)
#' 
#' @export
#' @import Seurat
#' 
transform_embedding <- function(
  obj,
  reduction = "umap",
  flip_horizontal = FALSE,
  flip_vertical = FALSE,
  rotate = 0
) {
  # Check if the reduction exists
  if (!reduction %in% names(obj@reductions)) {
    stop(paste("The reduction", reduction, "does not exist in the object"))
  }

  # Extract the embedding matrix
  embedding <- Embeddings(obj, reduction = reduction)
  if (is.null(embedding)) {
    stop(paste("The embedding", reduction, "does not exist in the object"))
  }

  # Flip the embedding
  if (flip_horizontal) {
    embedding[, 1] <- -embedding[, 1]
  }

  if (flip_vertical) {
    embedding[, 2] <- -embedding[, 2]
  }

  # Rotate the embedding
  if (rotate != 0) {
    theta <- rotate * (pi / 180)
    rotation_matrix <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2)
    embedding[, 1:2] <- embedding[, 1:2] %*% rotation_matrix
  }

  # Update the embedding matrix in the object
  obj[[reduction]]@cell.embeddings <- embedding

  # Return the object
  return(obj)
}