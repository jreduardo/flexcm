#-----------------------------------------------------------------------
#' @title COM-Poisson probability mass function
#' @param y a positive integer value.
#' @param mu the (approximate) mean parameter.
#' @param nu the dispersion parameter.
#' @param log logical. If \code{TRUE} it returns the logarithm of pmf.
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @export
#'
dcmp <- function(y, mu, nu, log = FALSE) {
  loglambda <- suppressWarnings(
    nu * log(mu + (nu - 1) / (2 * nu))
  )
  # Get the normalizing constants
  logz <- compute_logz(loglambda, mu, nu)
  # Compute the loglikelihood
  pmf <- y * loglambda - nu * lfactorial(y) - logz
  if (!log) pmf <- exp(pmf)
  return(pmf)
}

#-----------------------------------------------------------------------
#' @title Gamma-count probability mass function
#' @param y a positive integer value.
#' @param kappa the (asymptote) mean parameter.
#' @param alpha the dispersion parameter.
#' @param log logical. If \code{TRUE} it returns the logarithm of pmf.
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @export
#'
dgct <- function (y, kappa, alpha, log = FALSE) {
  alpha_kappa <- alpha * kappa
  alpha_y <- alpha * y
  pmf <- pgamma(1L, shape = alpha_y, rate = alpha_kappa) -
    pgamma(1L, shape = alpha_y + alpha, rate = alpha_kappa)
  if (log) pmf <- log(pmf)
  return(pmf)
}

#-----------------------------------------------------------------------
#' @title Discrete Weibull probability mass function
#' @param y a positive integer value.
#' @param q the "location" parameter (\code{0<q<1}).
#' @param rho the dispersion parameter.
#' @param log logical. If \code{TRUE} it returns the logarithm of pmf.
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @export
#'
ddwe <- function(y, q, rho, log = FALSE) {
  logq <- log(q)
  pmf <- exp(y^rho * logq) - exp((y + 1)^rho * logq)
  if (log) pmf <- log(pmf)
  return(pmf)
}

#-----------------------------------------------------------------------
#' @title Generalized Poisson probability mass function
#' @param y a positive integer value.
#' @param mu the mean parameter.
#' @param sigma the dispersion parameter.
#' @param log logical. If \code{TRUE} it returns the logarithm of pmf.
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @export
#'
dgpo <- function (y, mu, sigma, log = FALSE) {
  sigma_mu <- 1 + sigma * mu
  sigma_y <- 1 + sigma * y
  lsigma_mu <- suppressWarnings(log(sigma_mu))
  lsigma_y <- suppressWarnings(log(sigma_y))
  sy_sm <- suppressWarnings(sigma_y / sigma_mu)
  pmf <- y * (log(mu) - lsigma_mu) +
    (y - 1) * lsigma_y - mu *
    (sy_sm) - lfactorial(y)
  pmf[is.nan(pmf)] <- -Inf
  if (!log) pmf <- exp(pmf)
  return(pmf)
}

#-----------------------------------------------------------------------
#' @title Double Poisson probability mass function
#' @param y a positive integer value.
#' @param mu the (approximate) mean parameter.
#' @param phi the dispersion parameter.
#' @param log logical. If \code{TRUE} it returns the logarithm of pmf.
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @export
#'
ddpo <- function (y, mu, phi, log = FALSE) {
  # Get the normalizing constants
  logk <- compute_logk(mu, log(mu), phi, log(phi))
  ly <- log(y); ly[y == 0] <- 1
  pmf <- -0.5 * log(phi) - mu/phi - y + y * ly - lfactorial(y) +
    (y/phi) * (1 + log(mu) - ly) - logk
  if (!log) pmf <- exp(pmf)
  return(pmf)
}

#-----------------------------------------------------------------------
#' @title Poisson-Tweedie probability mass function
#' @param y a positive integer value.
#' @param mu the mean parameter.
#' @param omega the dispersion parameter.
#' @param power the power parameter.
#' @param npts number of points in the gaussian quadrature.
#' @param log logical. If \code{TRUE} it returns the logarithm of pmf.
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @importFrom stats dpois
#' @export
#'
dptw <- function(y, mu, omega, power, npts = 100L, log = FALSE) {
    pts <- statmod::gauss.quad(npts, kind = "laguerre")
    pjoint <- function(y, z) {
        dpois(y, lambda = z) *
            tweedie::dtweedie(z, mu = mu, phi = omega, power = power)
    }
    pmarginal <- function(y) {
        vapply(y, FUN = function(yi) {
            fpts <- pjoint(y = yi, z = pts$nodes)
            sum(pts$weights * fpts / exp(-pts$nodes))
        }, FUN.VALUE = numeric(1))
    }
    out <- vapply(y, FUN = function(yi) {
        if (yi == 0) {
            se <- sqrt(mu + omega * mu^power)
            yrange <- 1:min(1000, round(mu + 10 * se))
            probs <- pmarginal(yrange)
            out <- 1 - sum(probs)
        } else {
            pmarginal(yi)
        }
    }, FUN.VALUE = numeric(1))
    if (log) out <- log(out)
    return(out)
}
