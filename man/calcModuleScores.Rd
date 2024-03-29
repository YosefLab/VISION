% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-Module.R
\name{calcModuleScores}
\alias{calcModuleScores}
\title{Calculate module scores (signature scores but on the modules)}
\usage{
calcModuleScores(object, mod_norm_method = NULL, mod_gene_importance = TRUE)
}
\arguments{
\item{object}{the VISION object}

\item{mod_norm_method}{Method to apply to normalize the expression matrix
before calculating signature scores. Valid options are:
"znorm_columns" (default), "none", "znorm_rows", "znorm_rows_then_columns",
or "rank_norm_columns"}

\item{mod_gene_importance}{whether or not to rank each gene's contribution to
the overall signature score.  Default = TRUE.  This is used for inspecting
genes in a signature in the output report}
}
\value{
the VISION object, with the @ModScores and @ModGeneImportance slots populated
}
\description{
For each module-cell pair, compute a score that captures the level of
correspondence between the cell and the module.
}
