% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sichisq.R
\name{sichisq}
\alias{sichisq}
\alias{rsichisq}
\alias{dsichisq}
\title{Density and random sampling from the scaled inverse-chi-square distribution.}
\usage{
rsichisq(n, nu, tau)

dsichisq(x, nu, tau, log = FALSE)
}
\arguments{
\item{n}{Integer number of random draws to generate.}

\item{nu}{Shape parameter (positive scalar or vector).}

\item{tau}{Scale parameter (positive scalar or vector).}

\item{x}{Observations (positive scalar or vector).}

\item{log}{Logical; whether to calculate the density on the log scale.}
}
\description{
Density and random sampling from the scaled inverse-chi-square distribution.
}
\details{
Arguments \code{x}, \code{nu}, and \code{tau} are recycled as needed.
}
