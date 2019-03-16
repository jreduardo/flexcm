#-----------------------------------------------------------------------
#' @title Summary of the flexible count models (individual t-tests)
#' @param object an object of class \code{flexcm}, a result of call
#'   \code{\link{flexcm}(...)}.
#' @param correlation logical; if \code{TRUE}, the correlation matrix of
#'   the estimated parameters is returned and printed. Default is
#'   \code{FALSE}.
#' @param ... Currently not used.
#' @return an object of class \code{"summary.flexcm"}, a list with
#'   components.
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @importFrom stats cov2cor pnorm printCoefmat
#' @export
#'
summary.flexcm <- function(object, correlation = FALSE, ...) {
  if (is.null(object$vcov))
    stop(paste("Refit the model with `hessian=TRUE` to compute",
               "the standard errors."))
  power <- object$power
  powse <- attr(power, "std_error")
  #------------------------------------------
  estimates <- object$coefficients
  stderror  <- sqrt(diag(object$vcov))
  zvalue    <- estimates/stderror
  pvalue    <- 2 * pnorm(-abs(zvalue))
  ctable    <- cbind("Estimate"    = estimates,
                     "Std. Error"  = stderror,
                     "Z value"     = zvalue,
                     "Pr(>|z|)"    = pvalue)
  if (is.numeric(powse)) {
    pline <- c(power, powse, power/powse,
               2 * pnorm(-abs(power/powse)))
    ctable <- rbind("power" = pline, ctable)
  }
  #------------------------------------------
  if (correlation) {
    correlation <- cov2cor(object$vcov)
  }
  #--------------------------------------------
  out <- list(model       = object$model,
              coeftable   = ctable,
              loglik      = object$loglik,
              df.residual = object$df.residual,
              correlation = correlation,
              power       = object$power,
              call        = object$call)
  class(out) <- "summary.flexcm"
  return(out)
}

#-----------------------------------------------------------------------
# Print method for summary COM-Poisson models
#' @rdname flexcm-methods
#'
print.summary.flexcm <- function(x,
                                 digits = max(3L, getOption("digits") - 3L),
                                 ...) {
  name <- `_get_model_name`(x$model)
  index <- is.numeric(attr(x$power, "std_error"))
  cat(sprintf("\nIndividual Wald-tests for %s regression models",
              name), sep = "")
  cat("\nCall:  ",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  start_mean <- ifelse(index, 3, 2)
  cat("Mean coefficients:", "\n", sep = "")
  printCoefmat(x$coeftable[start_mean:nrow(x$coeftable), ,
                           drop = FALSE],
               digits = digits,
               has.Pvalue = TRUE)
  cat("\n")
  start_disp <- ifelse(index, 2, 1)
  cat("Dispersion coefficient:", "\n", sep = "")
  printCoefmat(x$coeftable[start_disp, , drop = FALSE],
               digits = digits,
               has.Pvalue = TRUE)
  cat("\n")
  if (name == "Poisson-Tweedie") {
    if (index) {
      cat("Power coefficient:", "\n", sep = "")
      printCoefmat(x$coeftable[1, , drop = FALSE],
                   digits = digits,
                   has.Pvalue = TRUE)
      cat("\n")
    } else {
      cat("Power coefficient (fixed): power = ", x$power,
          "\n\n", sep = "")
    }
  }
  if (is.matrix(x$correlation)) {
    cat("Correlation of coefficients:", "\n", sep = "")
    corr <- x$correlation
    corr <- format(round(corr, 2L), nsmall = 2L,
                   digits = digits)
    corr[!lower.tri(corr)] <- ""
    print(corr[-1, -ncol(corr), drop = FALSE], quote = FALSE)
    cat("\n")
  }
  cat("Residual degrees of freedom: ", x$df.residual,
      "\n", sep = "")
  cat("Minus twice the log-likelihood: ", -2 * x$loglik,
      "\n", sep = "")
  invisible(x)
}

#-----------------------------------------------------------------------
#' @name lrtests-flexcm
#' @title Likelihood ratio tests for nested flexible count models
#' @param object an object of class \code{flexcm}, a result of call
#'   \code{\link{flexcm}(...)}.
#' @param heading logical; if \code{TRUE}, the model formulas are
#'   printed as heading. Default is \code{TRUE}.
#' @param ... Currently not used.
#' @return an object of class \code{"anova"}, a table with compenents
#'   for LR test.
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @importFrom stats pchisq
#' @export
#'
anova.flexcm <- function(object, ..., heading) {
  dots <- list(...)
  isflexcm <- class(object) == "flexcm"
  #------------------------------------------
  # Organize the list of models
  if (isflexcm & length(dots) == 0)
    stop("Currently, this method only compare two (or more) models.")
  if (!isflexcm & length(dots) == 0) mlist <- object
  if (!isflexcm & length(dots) > 0)  mlist <- c(object, dots)
  if (isflexcm & length(dots)  > 0)  mlist <- c(list(object), dots)
  if (missing(heading)) heading <- is.null(names(mlist))
  #------------------------------------------
  # Handle the class of the models
  inherits_flexcm <- vapply(mlist, inherits, what = "flexcm", 0L)
  if (any(!inherits_flexcm))
    stop("Not all objects are of class \"flexcm\".")
  mnames <- vapply(mlist, function(m) m$model, "")
  if (!all(mnames[1L] == mnames)) {
    stop("The models were fitted under different likelihoods.")
  }
  obs <- vapply(mlist, function(x) x$nobs, 0L)
  if (diff(range(obs)) != 0) {
    warning("The models were fitted for different sample sizes.")
  }
  forms <- vapply(mlist, function(x) deparse(x$formula), "")
  #--------------------------------------------
  # Compute stats
  rds <- vapply(mlist, function(x) x$df.residual, 0L)
  lls <- vapply(mlist, function(x) x$loglik, 0)
  aic <- -2 * lls + 2 * (obs - rds)
  bic <- -2 * lls + log(obs) * (obs - rds)
  lrs <- c(NA, abs(diff(-2 * lls)))
  dfs <- c(NA, abs(diff(rds)))
  pvs <- pchisq(q = lrs, df = dfs, lower.tail = FALSE)
  nam <- if (is.null(names(mlist))) {
    sprintf("Model %i", seq_along(mlist))
  } else names(mlist)
  tab <- data.frame("Model"      = nam,
                    "AIC"        = aic,
                    "BIC"        = bic,
                    "Resid.df"   = rds,
                    "loglik"     = lls,
                    "Chisq.df"   = dfs,
                    "Chisq"      = lrs,
                    "Pr(>Chisq)" = pvs,
                    check.names = FALSE,
                    stringsAsFactors = FALSE)
  rownames(tab) <- NULL
  #--------------------------------------------
  # Heading
  if (heading) {
    diffpower <- vapply(mlist, function(x)
      typeof(attr(x$power, "std_error")), "")
    if (!all(diffpower[1] == diffpower)) {
      powers <- vapply(mlist, function(x) x$power, 0)
      heading <- paste(sprintf("  %s (power=%.2f): %s",
                               format(nam),
                               powers,
                               forms),
                       collapse = "\n")
    } else {
      heading <- paste(sprintf("  %s: %s",
                               format(nam),
                               forms),
                       collapse = "\n")
    }
  } else {
    heading <- NULL
  }
  name <- `_get_model_name`(mnames[1L])
  attr(tab, "name") <- name
  attr(tab, "heading") <- heading
  class(tab) <- c("anova.flexcm", "data.frame")
  return(tab)
}

#' @rdname lrtests-flexcm
#' @export
anova.list <- function(object, ..., heading) {
  out <- anova.flexcm(object, ..., heading = heading)
  return(out)
}

#' @export
print.anova.flexcm <- function(x,
                               digits = max(getOption("digits") - 2L, 3L),
                               signif.stars = getOption("show.signif.stars"),
                               ...) {
  name <- attr(x, "name")
  heading <- attr(x, "heading")
  cat(sprintf("\nLikelihood ratio test for %s", name),
      "regression models", "\n\n")
  if (!is.null(heading)) {
    cat(heading, "\n\n")
  }
  xx <- x
  rownames(xx) <- x$Model
  printCoefmat(xx[, -1],
               cs.ind = NA,
               tst.ind = c(3L, 5L),
               zap.ind = 7L,
               na.print = "",
               digits = digits,
               signif.stars = signif.stars)
  invisible(x)
}

#' @export
#' @rdname lrtests-flexcm
lrtest <- function (object, ..., heading) {
  UseMethod("anova", object)
}

#-----------------------------------------------------------------------
#' @title Likelihood ratio test for equidispersion
#' @param object an object (or a list of objects) of class
#'   \code{"flexcm"} with a model nested to Poisson (i.e. Com-Poisson,
#'   Gamma-count, Generalized-Poisson, Double Poisson and
#'   Poissson-Tweedie).
#' @param ... possibly, others fitted models (objects of class
#'   \code{"flexcm"}).
#' @param heading logical; if \code{TRUE}, the model formulas are
#'   printed as heading.
#' @details Note that, when more than one model is passed the models
#'   \strong{are not compared among themselves}. Here, each model is
#'   compared to Poisson and them organized in a table.
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @export
#'
equitest <- function (object, ..., heading) {
  UseMethod("equitest", object)
}

#' @export
equitest.flexcm <- function(object, ..., heading) {
  dots <- list(...)
  isflexcm <- class(object) == "flexcm"
  #------------------------------------------
  # Organize the list of models
  if (!isflexcm &  is.null(dots)) mlist <- object
  if ( isflexcm &  is.null(dots)) mlist <- list(object)
  if (!isflexcm & !is.null(dots)) mlist <- c(object, dots)
  if ( isflexcm & !is.null(dots)) mlist <- c(list(object), dots)
  if (any(!vapply(mlist, inherits, what = "flexcm", 0L)))
    stop("Not all objects are of class \"flexcm\".")
  if (missing(heading))
    heading <- ifelse(length(mlist) == 1,
                      FALSE,
                      is.null(names(mlist)))
  fam <- vapply(mlist, function(x) `_get_model_name`(x$model), "")
  if ("Discrete Weibull" %in% fam) {
    stop("Discrete Weibull models are not nested to Poisson.")
  }
  #--------------------------------------------
  lpo <- vapply(mlist, function(x) x$poissonfit$loglik, 0)
  npo <- vapply(mlist, function(x) length(x$poissonfit$coef), 0L)
  rds <- vapply(mlist, function(x) x$df.residual, 0L)
  lls <- vapply(mlist, function(x) x$loglik, 0)
  obs <- vapply(mlist, function(x) x$nobs, 0L)
  lrs <- 2 * (lls - lpo)
  dfs <- obs - rds - npo
  pvs <- pchisq(q = lrs, df = dfs, lower.tail = FALSE)
  nam <- if (is.null(names(mlist))) {
    sprintf("Model %i", seq_along(mlist))
  } else names(mlist)
  tab <- data.frame("Family"     = fam,
                    "Model"      = nam,
                    "Resid.df"   = rds,
                    "Loglik"     = lls,
                    "LRT_stat"   = lrs,
                    "LRT_df"     = dfs,
                    "Pr(>LRT_stat)" = pvs,
                    check.names = FALSE,
                    stringsAsFactors = FALSE)
  rownames(tab) <- NULL
  #--------------------------------------------
  if (heading) {
    calls <- vapply(mlist, function(x) deparse(x$call), "")
    cnames <- paste0("  ", format(nam), ": ")
    calls <- paste0(cnames, calls)
    calls <- paste(calls, collapse = "\n")
  } else {
    calls <- NULL
  }
  attr(tab, "heading") <- calls
  class(tab) <- c("equitest.flexcm", "data.frame")
  return(tab)
}

#' @export
equitest.list <- function(object, ..., heading) {
  out <- equitest.flexcm(object, ..., heading = heading)
  return(out)
}

#' @export
print.equitest.flexcm <- function(x,
                                  digits = max(getOption("digits") - 2L, 3L),
                                  signif.stars = getOption("show.signif.stars"),
                                  ...) {
  xx <- x
  heading <- attr(x, "heading")
  cat("\nLikelihood ratio test for equidispersion", "\n\n")
  if (!is.null(heading)) {
    if (nrow(x) == 1L) {
      heading <- gsub("Model 1:", "", heading)
    }
    cat(heading, "\n\n")
  }
  if (nrow(x) > 1L) {
    rownames(xx) <- paste(format(x$Family), format(x$Model))
  } else {
    rownames(xx) <- format(x$Family)
  }
  printCoefmat(xx[, -(1:2)],
               cs.ind = NA,
               tst.ind = c(2L, 3L),
               zap.ind = 5L,
               na.print = "",
               digits = digits,
               signif.stars = signif.stars,
               ...)
  invisible(x)
}
