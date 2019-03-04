#' @useDynLib flexcm, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#-----------------------------------------------------------------------
#' @title Negative log-likelihood function from COM-Poisson model
#' @description TODO: write the likelihood function.
#' @param params vector of model parameters.
#' @param X Design matrix related to the (approximate) mean parameter
#'   \eqn{\mu = \exp(X \beta)}.
#' @param y Vector of observed count data.
#' @return The computed negative log-likelihood function.
#' @author Eduardo Jr <edujrrib@gmail.com>
#'
llcmp <- function(params, X, y) {
  nu <- exp(params[1])
  beta <- params[-1]
  # Map the CMP parameters
  eta <- X %*% beta
  loglambda <- suppressWarnings(
    nu * log(exp(eta) + (nu - 1) / (2 * nu))
  )
  # Get the normalizing constants
  logz <- compute_logz(loglambda, nu)
  # Compute the loglikelihood
  ll <- sum(y * loglambda - nu * lfactorial(y) - logz)
  return(-ll)
}

#-----------------------------------------------------------------------
#' @title Negative log-likelihood function from Gamma-count model
#' @description TODO: write the likelihood function.
#' @param params vector of model parameters.
#' @param X Design matrix related to the (asymptotic) mean parameter
#'   \eqn{\kappa = \exp(X \beta)}.
#' @param y Vector of observed count data.
#' @return The computed negative log-likelihood function.
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @importFrom stats pgamma
#'
llgct <- function (params, X, y) {
  beta <- params[-1]
  alpha <- exp(params[1])
  # Map the GCT parameters
  eta <- X %*% beta
  kappa <- exp(eta)
  alpha_kappa <- alpha * kappa
  alpha_y <- alpha * y
  # Compute the loglikelihood
  ll <- sum(log(
    pgamma(1L, shape = alpha_y, rate = alpha_kappa) -
      pgamma(1L, shape = alpha_y + alpha, rate = alpha_kappa)
  ))
  return(-ll)
}

#-----------------------------------------------------------------------
#' @title Negative log-likelihood function from Discrete Weibull model
#' @description TODO: write the likelihood function.
#' @param params vector of model parameters.
#' @param X Design matrix related to the ("location") parameter
#'   \eqn{\log(-\log(q)) = \exp(X \beta)}.
#' @param y Vector of observed count data.
#' @return The computed negative log-likelihood function.
#' @author Eduardo Jr <edujrrib@gmail.com>
#'
lldwe <- function(params, X, y) {
  beta <- params[-1]
  rho <- exp(params[1])
  # Map the DWe parameters
  Xb <- X %*% beta
  logq <- -exp(Xb)
  # Compute the loglikelihood
  ll <- sum(log(exp(y^rho * logq) -
                  exp((y+1)^rho * logq)))
  return(-ll)
}

#-----------------------------------------------------------------------
#' @title Negative log-likelihood function from Generalized Poisson
#'   model
#' @description TODO: write the likelihood function.
#' @param params vector of model parameters.
#' @param X Design matrix related to the mean parameter \eqn{\mu =
#'   \exp(X \beta)}.
#' @param y Vector of observed count data.
#' @return The computed negative log-likelihood function.
#' @author Eduardo Jr <edujrrib@gmail.com>
#'
llgpo <- function (params, X, y) {
  beta <- params[-1]
  sigma <- params[1]
  # Map the GPo parameters
  xb <- X %*% beta
  mu <- exp(X %*% beta)
  sigma_mu <- 1 + sigma * mu
  sigma_y <- 1 + sigma * y
  # Suppress warnings
  lsigma_mu <- suppressWarnings(log(sigma_mu))
  lsigma_y <- suppressWarnings(log(sigma_y))
  sy_sm <- suppressWarnings(sigma_y / sigma_mu)
  # Compute the loglikelihood
  ll <- sum(y * (xb - lsigma_mu) +
              (y - 1) * lsigma_y -
                mu * (sy_sm) -
                lfactorial(y))
  return(-ll)
}

#-----------------------------------------------------------------------
#' @title Negative log-likelihood function from Double Poisson model
#'
#' @description TODO: write the likelihood function.
#' @param params vector of model parameters.
#' @param X Design matrix related to the mean parameter \eqn{\mu =
#'   \exp(X \beta)}.
#' @param y Vector of observed count data.
#' @return The computed negative log-likelihood function.
#' @author Eduardo Jr <edujrrib@gmail.com>
#'
lldpo <- function (params, X, y) {
  lphi <- params[1]
  beta <- params[-1]
  # Map the DPo parameters
  phi <- exp(lphi)
  eta <- X %*% beta
  mu <- exp(eta)
  # Get the normalizing constants
  logk <- compute_logk(mu, eta, phi, lphi)
  # Compute the loglikelihood
  ly <- log(y); ly[y == 0] <- 1
  ll <- sum(-0.5 * lphi - mu/phi - y + y*ly - lfactorial(y) +
              (y/phi) * (1 + eta - ly) - logk)
  return(-ll)
}
