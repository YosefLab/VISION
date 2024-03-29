% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utilities.R
\name{read_10x_h5_v3}
\alias{read_10x_h5_v3}
\title{Read 10x HDF5 Output - CellRanger 3.0}
\usage{
read_10x_h5_v3(h5_file, ensToSymbol = TRUE)
}
\arguments{
\item{h5_file}{path to h5 file}

\item{ensToSymbol}{bool denoting whether or not to perform label conversion}
}
\value{
a list with two items
        Expression: sparse count matrix with gene expression counts
        Antibody: sparse count matrix with antibody capture counts
}
\description{
Loads 10x output counts and converts expression to gene symbols
}
\details{
This version uses the h5 file produced by "cellranger count"

This file is typically in a folder that looks like:

    \code{outs/filtered_feature_bc_matrices_h5.h5}
}
