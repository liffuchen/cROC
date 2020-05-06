#' Extracts auc (se) and ROC curve estimates
#'
#' This function extracts auc (se) and ROC curve estimates from PHroc and PrevHroc model fits
#'
#' The proportional hazards ROC (PHroc) model assumes proportional hazards between healthy and diseased
#' populations: \eqn{h(y|d=1,x)/h(y|d=0,x) = \theta} and produces ROC curve \eqn{TPR = FPR^\theta} and AUC \eqn{1/(1+\theta)}.
#' With covariate, PHroc model takes the form of \eqn{h(y) = h_0(y) exp(b1 d + b2 x + b3 d*x)} where y, d, and x are score, status,
#' and covariate, respectively, implying \eqn{h(y|d=1,x)/h(y|d=0,x) = exp(b1 + b3 x)}. For binary x, \eqn{\theta = exp(b1)} and \eqn{exp(b1+b3)}
#' at \eqn{x = 0} and \eqn{x = 1}, respectively. The proportional reversed hazards ROC (PrevHroc) assumes proportional reversed hazards
#' instead, produces ROC curve \eqn{TPR = 1-(1-FPR)^{1/\theta}}, and shares the same AUC \eqn{1/(1+\theta)}.
#' With \eqn{\lambda(y)} as the reversed hazard function, PrevHroc model takes the form of \eqn{\lambda(y) = \lambda_0(y) exp(b1 d + b2 x + b3 d*x)},
#' implying \eqn{\lambda(y|d=1,x)/\lambda(y|d=0,x) = exp(b1 + b3 x)}. Same as in PHroc, \eqn{\theta = exp(b1)} and \eqn{exp(b1+b3)}
#' at \eqn{x = 0} and \eqn{x = 1}.
#'
#' Delta method is used to derive s.e. of auc.
#'
#' @param x1 a value of covariate x1 at which to estimate auc and ROC curve. Can be absent.
#' @param x2 a value of covariate x2 at which to estimate auc and ROC curve. Can be absent.
#' @param fit a coxph/croc object of model fit.
#' @param ROC logical, if TRUE, ROC estimates are produced in addition to AUC. Default is FALSE.
#' @param model a character string specifying the model for the conave ROC curve, one of PH and PrevH. Default PH.
#' @param ngrid number of equal-spaced FPR values on unit interval to compute TRP. Dafault 100.
#' @return either a vector of auc estimates or (if ROC = TRUE) a list of a vector of auc estimates and dataframe of FPR and TPR.
#' @author Zhen Chen (zhen.chen@@nih.gov)
#' @examples
#' auc(fit,"PH")
#' auc(x1=1,x2=.6,fit,"PrevH",ROC=TRUE)
#' @references
#' #' @export
auc <- function(x1 = 20, x2 = NULL, fit, model = "PH", ROC = FALSE, ngrid = 100) {
  k <- length(fit$coef)
  myformula <- switch((k+1)/2,"1/(1+exp(b1))","1/(1+exp(b1+x1*b3))","1/(1+exp(b1+x1*b4+x2*b5))","1/(1+exp(b1+x1*b4+x2*b5+x1*x2*b7))")
  auc <- car::deltaMethod(fit,myformula,parameterNames=paste("b",1:length(fit$coef),sep=""))
  row.names(auc) <- "auc"
  if(ROC) {
    if (model == "PH") {
       return(list(auc=auc,roc=data.frame(FPR=seq(0,1,l=100),TPR=seq(0,1,l=ngrid)**(1/auc$Est-1))))
    } else if (model == "PrevH") {
       return(list(auc=auc,roc=data.frame(FPR=seq(0,1,l=100),TPR=1-(1-seq(0,1,l=ngrid))**(1/(1/auc$Est-1)))))
    }
  } else {
    return(auc=auc)
  }
}
