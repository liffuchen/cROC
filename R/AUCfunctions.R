#' This function extracts auc (se) and ROC curve estimates from PHroc and PrHroc model fits
#'
#' xxx
#'
#' @param x1 a value of covariate x1 at which to estimate auc and ROC curve. Can be absent.
#' @param x2 a value of covariate x2 at which to estimate auc and ROC curve. Can be absent.
#' @param fit a coxph/croc object of model fit.
#' @param ROC logical, if TRUE, ROC estimates are produced in addition to AUC. Default is FALSE.
#' @param model a character string specifying the model for the conave ROC curve, one of PH and PrevH.
#' @param ngrid number of equal-spaced FPR values on unit interval to compute TRP.
#' @return either a vector of auc estimates or (if ROC = TRUE) a list of a vector of auc estimates and dataframe of FPR and TPR.
#' @export
auc <- function(x1=20,x2=NULL,fit,model = c("PH", "PrevH"),ROC=FALSE,ngrid=100) {
  model <- match.arg(model)
  k <- length(fit$coef)
  myformula <- switch((k+1)/2,"1/(1+exp(b1))","1/(1+exp(b1+x1*b3))","1/(1+exp(b1+x1*b4+x2*b5))","1/(1+exp(b1+x1*b4+x2*b5+x1*x2*b7))")
  auc <- car::deltaMethod(fit,myformula,parameterNames=paste("b",1:length(fit$coef),sep=""))
  if(ROC) {
    if (model == "PH") {
       return(list(auc=auc,roc=data.frame(FPR=seq(0,1,l=100),TPR=seq(0,1,l=ngrid)**(1/auc$Est-1))))
    } else if (model == "PrevH") {
       return(list(auc=auc,roc=data.frame(FPR=seq(0,1,l=100),TPR=1-(1-seq(0,1,l=ngrid))**(1/(1/auc$Est-1)))))
    }
  } else {
    return(auc)
  }
}
