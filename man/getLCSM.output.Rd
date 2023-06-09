% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LCSM.output_helper.R
\name{getLCSM.output}
\alias{getLCSM.output}
\title{Extract Point Estimates And Standard Errors of Latent Change Score Model with Time-invariant Covariates (If Any)}
\usage{
getLCSM.output(model, curveFun, growth_TIC, names)
}
\arguments{
\item{model}{An object representing a fitted latent change score model.}

\item{curveFun}{A string specifying the functional form of the growth curve. Supported options for latent change score
models include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"} (or \code{"EXP"}), \code{"Jenss-Bayley"}
(or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the value passed from \code{getLCSM()}.}

\item{growth_TIC}{A string or character vector specifying the column name(s) of time-invariant covariate(s) contributing to the
variability of growth factors if any. It takes the value passed from \code{getLCSM()}.}

\item{names}{A character vector specifying parameter names. It takes the value passed from \code{getLCSM()}.}
}
\value{
A data frame containing the point estimates and standard errors for parameters of a latent change  score
model with time-invariant covariates (if any).
}
\description{
This function computes and returns a data frame containing point estimates and standard errors for the parameters
of a latent change score model with time-invariant covariates (if any)
}
\keyword{internal}
