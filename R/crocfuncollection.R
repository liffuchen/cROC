#' Fit Proportional Reversed Hazards (PrevH) model
#'
#' This function fits the PrevH model given outcome y and covariates x
#'
#' @param x a matix of predictors. This should not include an intercept and already contains interactions and so on.
#' @param y a numeric vector containing the outcome (no censoring for now).
#' @return an object of class coxph representing the fit.
#' @export
prevh.fit <- function(x, y, ...)
{
  est <- survival::coxph(Surv((max(y)+y[1])-y)~x)
  return(est)
}

#' Fit Proportional Hazards (PH) model
#'
#' This function fits the PH model given outcome y and covariates x
#'
#' @param x a matix of predictors. This should not include an intercept and already contains interactions and so on.
#' @param y a numeric vector containing the outcome (no censoring for now).
#' @return an object of class coxph representing the fit.
#' @export
ph.fit <- function(x, y, ...)
{
  est <- survival::coxph(Surv(y)~x)
  return(est)
}

#' @export
croc <- function(x, ...) UseMethod("croc")

#' Estimate concave ROC curves and associated summary measures
#'
#' This function estimates concave curves and various summary measures.
#'
#' @param x a a matix of predictors. This should not include an intercept and already contains interactions and so on.
#' @param y a numeric vector containing the outcome (no censoring for now).
#' @param model a character string specifying the model for the conave ROC curve, one of PH, PrevH, PBN, Emprical.
#' @return an object of class coxph representing the fit.
#' @export
croc.default <- function(x,y,model, ...) {
  if (model %in% "PH") {
    fitter <- get("ph.fit")
  }
  else if (model %in% "PrevH") {
    fitter <- get("prevh.fit")
    # else if (model %in% "PBN") {
    #   fitter <- get("pbn.fit")
  }
  est <- fitter(x=x,y=y)
  return(est)
}

#' The formula wrapper for croc.default
#'
#' This function provides the formula wrapper for croc.default function.
#'
#' @param formula a formula object, with the response on the left of a ~ operator, and the terms on the right.
#' @param data a data.frame in which to interpret the variables named in the formula, or in the subset and the weights argument.
#' @param model a character string specifying the model for the conave ROC curve, one of PH, PrevH, PBN, Emprical.
#' @param subset an expression indicating which subset of the rows of data should be used in the fit. All observations are included by default.
#' @param weights a vector of case weights.
#' @return an object of class coxph and croc representing the fit.
#' @export
croc.formula <- function(formula, data, model=c("PH","PrevH"), subset, weights, na.action, ...)
{
  model <- match.arg(model)
  Call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  nf <- names(mf)
  m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (nrow(mf) == 0) stop("No observations in data")
  Terms <- terms(mf)

  Y <- model.extract(mf, "response")
  contrast.arg <- NULL
  attr(Terms, "intercept") <- 1
  X <- model.matrix(Terms, mf, contrasts = contrast.arg)
  Xatt <- attributes(X)
  adrop <- c(0, untangle.specials(Terms, "strata")$terms)
  xdrop <- Xatt$assign %in% adrop
  X <- X[, !xdrop, drop = FALSE]

  fit <- croc.default(X,Y,model)

  fit$call <- Call
  class(fit) <- c("croc","coxph")
  return(fit)
}

#' Print method for class "croc"
#'
#' This function provides the print method for class "croc".
#'
#' @param x an object of class "croc".
#' @param digits the number of significant digits to use when printing.
#' @param signif.stars logical. If TRUE, ‘significance stars’ are printed for each coefficient.
#' @param getOption ????.
#' @return A list of statistics????.
#' @export
print.croc <- function(x, digits=max(1L, getOption("digits") - 3L), signif.stars=FALSE, ...) {
  if (!is.null(cl<- x$call)) {
    cat("Call:\n")
    dput(cl)
    cat("\n")
  }
  if (!is.null(x$fail)) {
    cat(" croc failed.", x$fail, "\n")
    return()
  }
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  coef <- x$coefficients
  names(coef) <- sub("^.","",names(coef))
  se <- sqrt(diag(x$var))
  if(is.null(coef) | is.null(se)) stop("Input is not valid")

  if (is.null(x$naive.var)) {
    tmp <- cbind(coef, exp(coef), se, coef/se, pchisq((coef/se)^2, 1, lower.tail=FALSE))
    dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)","se(coef)", "z", "p"))
  }
  else {
    nse <- sqrt(diag(x$naive.var))
    tmp <- cbind(coef, exp(coef), nse, se, coef/se, pchisq((coef/se)^2, 1, lower.tail=FALSE))
    dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)", "se(coef)", "robust se", "z", "p"))
  }
  printCoefmat(tmp, digits=digits, P.values=TRUE, has.Pvalue=TRUE,signif.stars = signif.stars, ...)

  logtest <- -2 * (x$loglik[1] - x$loglik[2])
  if (is.null(x$df)) df <- sum(!is.na(coef))
  else  df <- round(sum(x$df),2)
  cat("\n")
  cat("Likelihood ratio test=", format(round(logtest, 2)), "  on ",
      df, " df,", " p=",
      format.pval(pchisq(logtest, df, lower.tail=FALSE), digits=digits),
      "\n",  sep="")
  omit <- x$na.action
  cat("Total sample size:", x$n)
  if (length(omit)) cat("\   (", naprint(omit), ")\n", sep="")
  invisible(x)
}
