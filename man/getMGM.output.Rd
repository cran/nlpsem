% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MGM.output_helper.R
\name{getMGM.output}
\alias{getMGM.output}
\title{Extract Point Estimates And Standard Errors of Multivariate Latent Growth Curve Models Or Multivariate Latent Change
Score Models}
\usage{
getMGM.output(model, y_var, records, curveFun, y_model, names)
}
\arguments{
\item{model}{An object representing a fitted multivariate latent growth curve model or latent change score model.}

\item{y_var}{A vector of strings, with each element representing the prefix for column names corresponding to a
particular outcome variable at each study wave. It takes the value passed from \code{getMGM()}.}

\item{records}{A list of numeric vectors, with each vector specifying the indices of the observed study waves for
the corresponding outcome variable. It takes the value passed from \code{getMGM()}.}

\item{curveFun}{A string specifying the functional form of the growth curve. Supported options for \code{y_model =
"LGCM"} include: \code{"linear"} (or \code{"LIN"}), \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
(or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"bilinear spline"} (or \code{"BLS"}). Supported
options for \code{y_model = "LCSM"} include: \code{"quadratic"} (or \code{"QUAD"}), \code{"negative exponential"}
(or \code{"EXP"}), \code{"Jenss-Bayley"} (or \code{"JB"}), and \code{"nonparametric"} (or \code{"NonP"}). It takes the
value passed from \code{getMGM()}.}

\item{y_model}{A string specifying how to fit the longitudinal outcome. Supported values are \code{"LGCM"} and \code{"LCSM"}.
It takes the value passed from \code{getMGM()}.}

\item{names}{A character vector specifying parameter names. It takes the value passed from \code{getMGM()}.}
}
\value{
A data frame containing the point estimates and standard errors for parameters of a multivariate
latent growth curve model or a multivariate latent change score model.
}
\description{
This function computes and returns a data frame containing point estimates and standard errors for the parameters of
a multivariate latent growth curve model or a multivariate latent change score model.
}
\keyword{internal}
