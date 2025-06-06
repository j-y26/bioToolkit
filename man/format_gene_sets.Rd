% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/module_scoring_rnaseq.R
\name{format_gene_sets}
\alias{format_gene_sets}
\title{Manipulating geneset formats}
\usage{
format_gene_sets(gsc, to = c("GeneSetCollection", "list", "gmt"))
}
\arguments{
\item{gsc}{A GeneSetCollection object or a list of gene sets or a path to a GMT file.}

\item{to}{A character string indicating the format to convert to. Available
options are "GeneSetCollection", "list", and "gmt".}
}
\value{
A GeneSetCollection object, a list of gene sets, or a GMT file.
}
\description{
This function is used to manipulate gene set formats
             such as GMT files and GeneSetCollection objects.
}
