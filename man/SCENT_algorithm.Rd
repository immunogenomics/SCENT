% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCENTfunctions.R
\name{SCENT_algorithm}
\alias{SCENT_algorithm}
\title{SCENT Algorithm: Poisson Regression with Empirical P-values through Bootstrapping.}
\usage{
SCENT_algorithm(object, celltype, ncores, regr = "poisson", bin = TRUE)
}
\arguments{
\item{object}{SCENT object}

\item{celltype}{character. User specified cell type defined in celltypes column of meta.data}

\item{ncores}{numeric. Number of cores to use for Parallelization}

\item{regr}{character. Regression type: "poisson" or "negbin" for Poisson regression and Negative Binomial regression, respectively}

\item{bin}{logical. TRUE to binarize ATAC counts. FALSE to NOT binarize ATAC counts}
}
\value{
SCENT object with updated field SCENT.results
}
\description{
SCENT Algorithm: Poisson Regression with Empirical P-values through Bootstrapping.
}
