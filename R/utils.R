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
