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
  mu <- exp(eta)
  loglambda <- suppressWarnings(
    nu * log(exp(eta) + (nu - 1) / (2 * nu))
  )
  # Get the normalizing constants
  logz <- compute_logz(loglambda, mu, nu)
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

#-----------------------------------------------------------------------
#' @title Maximize log-likelihood functions
#' @param llfun Function thar returns the negative of the log-likelihood
#'   given \code{params}, \code{X} and \code{y} arguments for some count
#'   model.
#' @param X Design matrix.
#' @param y Vector of observed count data.
#' @param start Initial parameters
#' @param method Argument passed to \code{\link[stats]{optim}(...)}.
#' @param lower Argument passed to \code{\link[stats]{optim}(...)}.
#' @param upper Argument passed to \code{\link[stats]{optim}(...)}.
#' @param hessian Argument passed to \code{\link[stats]{optim}(...)}.
#' @param control Argument passed to \code{\link[stats]{optim}(...)}.
#' @return A list with estimated parameters and hessian matrix.
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @importFrom stats glm.fit poisson optim
#'
flexcm_fit <- function(llfun, X, y,
                       start   = NULL,
                       method  = c("BFGS",
                                   "Nelder-Mead",
                                   "CG",
                                   "L-BFGS-B",
                                   "SANN"),
                       lower   = -Inf,
                       upper   = Inf,
                       hessian = TRUE,
                       control = list()) {
  #-------------------------------------------
  # Initial values
  if (is.null(start)) {
    model <- glm.fit(x = X, y = y, family = poisson())
    start <- c("disp" = 0, model$coefficients)
  } else {
    if (is.null(names(start)))
      names(start)  <- c("disp", colnames(X))
  }
  #------------------------------------------
  # Maximization
  method <- match.arg(method)
  out <- optim(par = start,
               fn = llfun,
               method = method,
               lower = lower,
               upper = upper,
               hessian = hessian,
               control = control,
               X = X,
               y = y)
  return(out)
}

#-----------------------------------------------------------------------
#' @title Fitting flexible count model models
#' @param formula An object of class "\code{\link[stats]{formula}}"
#'   describe the model for mean.
#' @param data An optional data frame containing the variables in the
#'   model. If not found in data, the variables are taken from
#'   \code{environment(formula)}.
#' @param model An character indicating the model used. Options are
#'   \code{"compoisson"}, \code{"gammacount"}, \code{"discreteweibull"},
#'   \code{"generalizedpoisson"}, \code{"doublepoisson"}, and
#'   \code{"poissontweedie"}.
#' @param ... Arguments to be used by \code{\link{flexcm_fit}}.
#' @return An object of class \code{cmpreg}.
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @importFrom stats model.frame model.matrix model.response setNames
#' @export
#'
flexcm <- function(formula, data, model, ...) {
  #--------------------------------------------
  if (missing(data))
    data <- environment(formula)
  #-------------------------------------------
  # Get matrices
  frame <- model.frame(formula, data)
  terms <- attr(frame, "terms")
  X <- model.matrix(terms, frame)
  y <- model.response(frame)
  #-------------------------------------------
  # Get the log-likelihood function
  llfun <- switch(model,
                  "compoisson"         = llcmp,
                  "gammacount"         = llgct,
                  "discreteweibull"    = lldwe,
                  "generalizedpoisson" = llgpo,
                  "doublepoisson"      = lldpo,
                  "poissontweedie"     = NA,
                  # Abbreviation
                  "cmp"                = llcmp,
                  "gct"                = llgct,
                  "dwe"                = lldwe,
                  "gpo"                = llgpo,
                  "dpo"                = lldpo,
                  "ptw"                = NA)
  dname <- switch(model,
                  "compoisson"         = "log(nu)",
                  "gammacount"         = "log(alpha)",
                  "discreteweibull"    = "log(rho)",
                  "generalizedpoisson" = "sigma",
                  "doublepoisson"      = "log(phi)",
                  "poissontweedie"     = "omega",
                  # Abbreviation
                  "cmp"                = "log(nu)",
                  "gct"                = "log(alpha)",
                  "dwe"                = "log(rho)",
                  "gpo"                = "sigma",
                  "dpo"                = "log(phi)",
                  "ptw"                = "omega")
  #--------------------------------------------
  # Optimize
  details <- flexcm_fit(llfun = llfun, X = X, y = y, ...)
  coefficients <- setNames(details$par, c(dname, colnames(X)))
  mean_coefficient <- coefficients[-1]
  disp_coefficient <- coefficients[ 1]
  vcov <- NULL
  if ("hessian" %in% names(details))
    vcov <- solve(details$hessian)
  #--------------------------------------------
  # Fitted values
  fitted <- exp(X %*% mean_coefficient)
  if (model %in% c("discreteweibull", "dwe"))
    fitted <- exp(-fitted)
  #--------------------------------------------
  # Output
  out <- list(call = match.call(),
              model = model,
              formula = formula,
              nobs = length(y),
              df.residual = length(y) - length(details$par),
              details = details,
              loglik = -details$value,
              vcov = vcov,
              coefficients = coefficients,
              mean_coefficient = mean_coefficient,
              disp_coefficient = disp_coefficient,
              fitted = c(fitted),
              data = list(X = X, y = y))
  class(out) <- "flexcm"
  return(out)
}
