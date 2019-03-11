#-----------------------------------------------------------------------
# Get the name of the models
`_get_model_name` <- function(model){
  if (model %in% c("compoisson", "cmp")) {
    "COM-Poisson"
  } else if ((model %in% c("gammacount", "gct"))) {
    "Gamma-count"
  } else if ((model %in% c("discreteweibull", "dwe"))) {
    "Discrete Weibull"
  } else if ((model %in% c("generalizedpoisson", "gpo"))) {
    "Generalized Poisson"
  } else if ((model %in% c("doublepoisson", "dpo"))) {
    "Double Poisson"
  } else if ((model %in% c("poissontweedie", "ptw"))) {
    "Poisson-Tweedie"
  } else {
    NULL
  }
}

#-----------------------------------------------------------------------
# Compute the loglikelihood for Poisson-Tweedie models
`_ptw_loglik` <- function(y, mu, omega, power, npts = 100L) {
  # Conditions
  if (!(omega >= 0 & (power >= 1 | power <= 0))) {
    ll <- NA
  } else {
    lli <- vapply(seq(mu), function(i)
      dptw(y = y[i], mu = mu[i], omega = omega,
           power = power, npts = npts, log = TRUE),
      numeric(1))
    ll <- sum(lli)
  }
  return(ll)
}

#-----------------------------------------------------------------------
# Compute the expectation of Gamma-count
`_compute_mean_gct` <- function(eta, dispersion, tol = 1e-5) {
  kappa <- exp(eta)
  alpha <- exp(dispersion)
  #------------------------------------------
  ymax <- ceiling(kappa + 5 * sqrt(kappa/alpha))
  prob <- dgct(ymax, kappa, alpha)
  if (any(prob > tol)) {
    for (i in which(prob > tol)) {
      for (j in 1:100) {
        ymax[i] <- ymax[i] + 5
        prob[i] <- dgct(ymax[i], kappa[i], alpha)
        if (prob[i] < tol) break
      }
    }
  }
  out <- mapply(ymax, kappa, SIMPLIFY = FALSE,
                FUN = function(ymaxi, kappai) {
                  range <- 1:ymaxi
                  # sum(range * dgct(range, kappai, alpha))
                  sum(pgamma(alpha * kappai, range * alpha))
                })
  unlist(out)
}

#-----------------------------------------------------------------------
# Compute the expectation of Discrete Weibull
`_compute_mean_dwe` <- function(eta, dispersion, tol = 1e-5) {
  q <- exp(-exp(eta))
  rho <- exp(dispersion)
  #--------------------------------------------
  lambda <- -log(q)
  scale = exp(-log(lambda)/rho)
  approx_me <- scale * gamma(1 + 1/rho)
  approx_va <- scale^2 * (gamma(1 + 2/rho) - (gamma(1 + 1/rho))^2)
  ymax <- ceiling(approx_me + 5 * sqrt(approx_va))
  prob <- ddwe(ymax, q, rho)
  if (any(prob > tol)) {
    for (i in which(prob > tol)) {
      for (j in 1:100) {
        ymax[i] <- ymax[i] + 5
        prob[i] <- ddwe(ymax[i], q[i], rho)
        if (prob[i] < tol) break
      }
    }
  }
  out <- mapply(ymax, q, SIMPLIFY = FALSE,
                FUN = function(ymaxi, qi) {
                  range <- 1:ymaxi
                  # sum(range * ddwe(range, qi, rho))
                  sum(qi^range^rho)
                })
  unlist(out)
}
