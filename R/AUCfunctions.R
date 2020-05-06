#' Extacts auc (s.e.) and ROC curve estimates from PHroc and PrHroc model fits
#'
#' This function extracts auc (se) and ROC curve estimates from PHroc and PrHroc model fits with <= 2 covariates.
#'
#' @param x1 a value of covariate x1 at which to estimate auc and ROC curve. Can be absent.
#' @param x2 a value of covariate x2 at which to estimate auc and ROC curve. Can be absent.
#' @param fit a coxph/croc object of model fit.
#' @param ROC logical, if TRUE, ROC estimates are produced in addition to AUC. Default is FALSE.
#' @param ngrid number of equal-spaced FPR values on unit interval to compute TRP.
#' @return either a vector of auc estimates or (if ROC = TRUE) a list of a vector of auc estimates and dataframe of FPR and TPR.
#' @export
auc <- function(x1=20,x2=NULL,fit,ROC=FALSE,ngrid=100) {
  k <- length(fit$coef)
  myformula <- switch((k+1)/2,"1/(1+exp(b1))","1/(1+exp(b1+x1*b3))","1/(1+exp(b1+x1*b4+x2*b5))","1/(1+exp(b1+x1*b4+x2*b5+x1*x2*b7))")
  auc <- car::deltaMethod(fit,myformula,parameterNames=paste("b",1:length(fit$coef),sep=""))
  if(ROC) {
    return(list(auc=auc,roc=data.frame(FPR=seq(0,1,l=100),TPR=seq(0,1,l=ngrid)**(1/auc$Est-1))))
  } else {
    return(auc)
  }
}
