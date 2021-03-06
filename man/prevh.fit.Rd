% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crocfuncollection.R
\name{prevh.fit}
\alias{prevh.fit}
\title{Fit Proportional Reversed Hazards (PrevH) model}
\usage{
prevh.fit(x, y, ...)
}
\arguments{
\item{x}{a matix of predictors. This should not include an intercept.}

\item{y}{a numeric vector containing the outcome (no censoring for now).}
}
\value{
an object of class coxph representing the fit.
}
\description{
This function fits the PrevH model given outcome y and covariates x
}
