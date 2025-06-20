% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/single_cell_visualize.R
\name{grouped_bar_plot}
\alias{grouped_bar_plot}
\title{Bar Plot for Grouped Data}
\usage{
grouped_bar_plot(
  df,
  group.by,
  split.by,
  frequency,
  mode = c("stack", "fill"),
  palette = "polychrome",
  title = NULL,
  show_totals = FALSE
)
}
\arguments{
\item{df}{A data frame with columns: `split.by`, `group.by`, and `frequency`.}

\item{group.by}{Column name used for the fill (segments within bars).}

\item{split.by}{Column name used for the x-axis (categories to split the bars).}

\item{frequency}{Column name with the values for bar heights (usually counts 
or proportions).}

\item{mode}{Either "stack" or "fill".}

\item{palette}{A valid palette name from "alphabet", "alphabet2", "glasbey", "polychrome", "stepped", and "parade".
Alternatively, a vector of colors can be provided with the same length as unique groups.}

\item{title}{Plot title.}

\item{show_totals}{Logical, whether to show totals on top of the bars.}
}
\value{
A ggplot2 object
}
\description{
Bar Plot for Grouped Data
}
