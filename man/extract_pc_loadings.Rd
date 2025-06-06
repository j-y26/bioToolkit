% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_rnaseq.R
\name{extract_pc_loadings}
\alias{extract_pc_loadings}
\title{Extract PCA loadings from PCA object}
\usage{
extract_pc_loadings(pca_obj, pcs = 1:3, top_n = 20)
}
\arguments{
\item{pca_obj}{A PCA object generated using the \code{prcomp} function.}

\item{pcs}{A vector of integers specifying the principal components to
extract the loadings. Default is to extract PC 1-3.}

\item{top_n}{The number of top genes to extract for each principal component.
Default is 20.}
}
\value{
A list, where each element is a PC, and each element is a named vector
        containing the top genes and their loadings.
}
\description{
This function extracts the PCA loadings from a PCA object
             generated using the \code{prcomp} function.
}
