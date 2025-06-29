% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/single_cell_atac.R
\name{peak_category_qc}
\alias{peak_category_qc}
\title{Peak category QC across one or more peak calling methods/samples}
\usage{
peak_category_qc(peak_list, peak_names = NULL, TxDb)
}
\arguments{
\item{peak_list}{A named list of GRanges objects, each representing peaks
identified by a different peak calling method or sample.}

\item{peak_names}{A character vector of names corresponding to the peaks in
`peak_list`. If not provided, names will be extracted from
`peak_list` names. Default is `NULL`.}

\item{TxDb}{A TxDb object containing gene annotations. This is used to
categorize peaks into promoter, intronic, intergenic, and
exonic regions.}
}
\value{
A list of csAnno objects, each containing the detailed annotation of
        each peak set. See `ChIPseeker` package for more details.
}
\description{
This function performs quality control on peak categories across
             one or more peak calling methods. It generates a bar plot
             showing the distribution of peaks across different categories.
}
