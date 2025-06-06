% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/single_cell_visualize.R
\name{add_dimension_arrows}
\alias{add_dimension_arrows}
\title{Add dimension arrows}
\usage{
add_dimension_arrows(
  plot,
  reduction = "umap",
  dims = NULL,
  arrow_length = 2,
  x_coord_adjust = 0.25,
  y_coord_adjust = x_coord_adjust,
  arrow_size = 0.25,
  arrow_width = 0.2,
  arrow_line_width = 3,
  arrow_color = "black",
  text_font_size = 12,
  text_dst = 0.3
)
}
\arguments{
\item{plot}{A ggplot object}

\item{reduction}{The reduction method used to generate the plot, default is
"umap", allowed values are "umap", "tsne", "pca", "spatial"}

\item{dims}{A vector of dimension names, overwrites the default dimension
names assumed from reduction method}

\item{arrow_length}{The size of the arrow, default is 2 (cm)}

\item{x_coord_adjust}{Horizontal adjustment of the arrow origin, default is 0.25 (cm)}

\item{y_coord_adjust}{Vertical adjustment of the arrow origin, default is the same as x_coord_adjust}

\item{arrow_size}{The size of the arrow head, default is 0.25 cm}

\item{arrow_width}{The width of the arrow, default is 0.2 cm}

\item{arrow_line_width}{The width of the arrow line, default is 3}

\item{arrow_color}{The color of the arrow and label, default is "black"}

\item{text_font_size}{The font size of the text, default is 12}

\item{text_dst}{The distance of the text from the arrow in cm, default is 0.3}
}
\value{
a ggplot object with the dimension arrows added
}
\description{
This function adds dimension arrows to a plot at the bottom
             left corner of the plot. The arrows represent the direction of
             the dimensions in the plot. Assumes that the plot is a ggplot
             object and does not have an existing axis.
}
