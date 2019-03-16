#-----------------------------------------------------------------------
#' @title Predict method for flexible count models
#' @details TODO
#' @param object a fitted object of class from \code{"flexcm"}.
#' @param newdata optionally, a data frame in which to look for
#'   variables with which to predict. If omitted, the fitted linear
#'   predictors are used.
#' @param newmatrix optionally, a list with named design matrices
#'   (\code{"X"}) used to predict. If omitted, the fitted linear
#'   predictors are used.
#' @param type the type of prediction required. The default
#'   \code{"link"} is on the scale of the linear predictors; the
#'   alternative \code{"response"} is on the scale of the response
#'   variable (for the models that the location parameter not is the
#'   mean, it will be computed the mean by \eqn{\sum x p(x; \theta)}).
#' @param interval character indicating the type of interval
#'   calculation. Currently the intervals can be \code{"none"} or
#'   \code{"confidence"}.
#' @param level level of confidence intervals.
#' @param augment_data logical indicating if \code{newdata} should be
#'   augmented with the predict values and, possibly, confidence
#'   intervals. Default is \code{FALSE}.
#' @param ... currently not used.
#' @return a tibble with fitted values (column code{fit}) and, possibly,
#'   confidence intervals (columns \code{lwr} and \code{upr}). If
#'   \code{augment_data}, the result will be the \code{newdata} with new
#'   columns \code{fit}, fitted values, and code{ste}, standard errors.
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @importFrom stats qnorm
#' @export
#'
predict.flexcm <- function(object,
                           newdata,
                           newmatrix = NULL,
                           type = c("link", "response"),
                           interval = c("none", "confidence"),
                           level = 0.95,
                           augment_data = FALSE,
                           ...) {
  mname <- `_get_model_name`(object$model)
  #------------------------------------------
  type <- match.arg(type)
  interval <- match.arg(interval)
  missingdata <- missing(newdata)
  missingmatx <- is.null(newmatrix)
  #------------------------------------------
  Vcov <- object$vcov
  Vbeta <- Vcov[-1, -1, drop = FALSE]
  Vdisp <- Vcov[ 1,  1, drop = FALSE]
  Vbedi <- Vcov[-1,  1, drop = FALSE]
  #------------------------------------------
  if (!missingdata & !missingmatx)
    stop("Use only 'newdata' or 'newmatrices'.")
  if ( missingdata &  missingmatx)
    X <- object$data$X
  if (!missingmatx &  missingdata)
    X <- newmatrix
  if ( missingmatx & !missingdata)
    X <- model.matrix(object$formula[-2], newdata)
  #------------------------------------------
  if (mname %in% "Poisson-Tweedie")
    Vcond <- Vbeta
  Vcond <- Vbeta - tcrossprod(Vbedi %*% solve(Vdisp), Vbedi)
  beta <- object$mean_coefficient
  out <- cbind(fit = c(X %*% beta))
  #--------------------------------------------
  if (interval == "confidence") {
    qn <- -qnorm((1 - level)/2)
    ste <- sqrt(diag(tcrossprod(X %*% Vcond, X)))
    out <- cbind(out,
                 out[, "fit"] - qn * ste,
                 out[, "fit"] + qn * ste)
    if (mname == "Discrete Weibull") out <- out[, c(1, 3, 2)]
    colnames(out) <- c("fit", "lwr", "upr")
  }
  #------------------------------------------
  if (type == "response") {
    func <- switch(mname,
                   "Gamma-count"      = `_compute_mean_gct`,
                   "Discrete Weibull" = `_compute_mean_dwe`,
                   function(eta, dispersion) exp(eta))
    out <- apply(out, 2, func,
                 dispersion = object$disp_coefficient)
  }
  if (augment_data & missingdata)
    warning("Needs 'newdata' for augment data.")
  if (augment_data & !missingdata)
    out <- cbind(newdata, out)
  out <- as.data.frame(out)
  class(out) <- c("tbl_df", "tbl", "data.frame")
  return(out)
}

#-----------------------------------------------------------------------
# Fitted method
#' @rdname flexcm-methods
#' @export
#'
fitted.flexcm <- function(object, ...) {
  object$fitted
}
