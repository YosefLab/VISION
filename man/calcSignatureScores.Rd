% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnalysisFunctions.R
\name{calcSignatureScores}
\alias{calcSignatureScores}
\title{calculate signature scores}
\usage{
calcSignatureScores(object, sig_norm_method = NULL, sig_gene_importance = TRUE)
}
\arguments{
\item{object}{the VISION object}

\item{sig_norm_method}{Method to apply to normalize the expression matrix
before calculating signature scores. Valid options are:
"znorm_columns" (default), "none", "znorm_rows", "znorm_rows_then_columns",
or "rank_norm_columns"}

\item{sig_gene_importance}{whether or not to rank each gene's contribution to
the overall signature score.  Default = TRUE.  This is used for inspecting
genes in a signature in the output report}
}
\value{
the VISION object, with the @SigScores and @SigGeneImportance slots populated
}
\description{
For each signature-cell pair, compute a score that captures the level of
correspondence between the cell and the signature.
}
