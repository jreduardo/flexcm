#' @name flexcm-methods
#' @title Method for \code{'flexcm'} objects
#' @param object an object of class \code{"flexcm"}.
#' @param x an object of class \code{"flexcm"}.
#' @param digits minimal number of _significant_ digits, see
#'   \code{\link[base]{print.default}}.
#' @param ... Currently not used.
#' @return .
#' @author Eduardo Jr <edujrrib@gmail.com>
#'
NULL

#-----------------------------------------------------------------------
# Print method
#' @rdname flexcm-methods
#' @export
#'
print.flexcm <- function(x,
                         digits = max(3L, getOption("digits") - 3L),
                         ...) {
  name <- `_get_model_name`(x$model)
  cat(sprintf("\n%s regression models", name), sep = "")
  cat("\nCall:  ",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Mean coefficients:", "\n", sep = "")
  print.default(format(x$mean_coefficient, digits = digits),
                print.gap = 2, quote = FALSE)
  cat("\n")
  cat(sprintf("Dispersion coefficient: %s = ",
              names(x$disp_coefficient)),
      format(x$disp_coefficient, digits = digits),
      "\n", sep = "")
  if (name == "Poisson-Tweedie") {
    aux <- ifelse(is.numeric(attr(x$power, "std_error")),
                  "estimated", "fixed")
    cat(sprintf("Power coefficient (%s): power = ", aux),
        format(x$power, digits = digits),
        "\n", sep = "")
  }
  cat("\n")
  cat("Residual degrees of freedom: ", x$df.residual,
      "\n", sep = "")
  cat("Minus twice the log-likelihood: ", -2 * x$loglik,
      "\n", sep = "")
  invisible(x)
}

#-----------------------------------------------------------------------
# Get the log-likelihood
#' @rdname flexcm-methods
#' @export
#'
logLik.flexcm <- function(object, ...) {
  if (!missing(...))
    warning("extra arguments discarded")
  ll <- object$loglik
  attr(ll, "df") <- object$nobs - object$df.residual
  attr(ll, "nobs") <- object$nobs
  class(ll) <- "logLik"
  return(ll)
}

#-----------------------------------------------------------------------
# Get the parameter estimates
#' @rdname flexcm-methods
#' @export
#'
coef.flexcm <- function(object, ...) {
  if (!missing(...))
    warning("extra arguments discarded")
  object$coefficients
}

#-----------------------------------------------------------------------
# Get the variance-covariance matrix
#' @rdname flexcm-methods
#' @export
#'
vcov.flexcm <- function(object, ...) {
  if (!missing(...))
    warning("extra arguments discarded")
  object$vcov
}

#-----------------------------------------------------------------------
# Get the design matrices
#' @rdname flexcm-methods
#' @export
#'
model.matrix.flexcm <- function(object, ...) {
  if (!missing(...))
    warning("extra arguments discarded")
  object$data$X
}
