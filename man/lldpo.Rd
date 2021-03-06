% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihoods.R
\name{lldpo}
\alias{lldpo}
\title{Negative log-likelihood function from Double Poisson model}
\usage{
lldpo(params, X, y)
}
\arguments{
\item{params}{vector of model parameters.}

\item{X}{Design matrix related to the mean parameter \eqn{\mu =
\exp(X \beta)}.}

\item{y}{Vector of observed count data.}
}
\value{
The computed negative log-likelihood function.
}
\description{
TODO: write the likelihood function.
}
\author{
Eduardo Jr <edujrrib@gmail.com>
}
