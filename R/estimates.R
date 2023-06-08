#---------------------------------------------
#' Variance-covariance matrix
#'
#' @aliases vcov.survstan
#' @description This function extracts and returns the variance-covariance matrix associated with the regression coefficients when the maximum likelihood estimation approach is used in the model fitting.
#' @export
#' @importFrom stats vcov
#' @param object an object of the class survstan.
#' @param all logical; if FALSE (default), only covariance matrix associated with regression coefficients is returned; if TRUE, the full covariance matrix is returned.
#' @param ... further arguments passed to or from other methods.
#' @return  the variance-covariance matrix associated with the parameters estimators.
#' @examples
#' \donttest{
#' library(survstan)
#' fit <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
#' vcov(fit)
#' }
#'
vcov.survstan <- function(object, all = FALSE, ...){
  if(object$survreg == "yp"){
    p <- 2*object$p
  }else{
    p <- object$p
  }

  labels <- names(object$estimates)
  if(!isTRUE(all) & p > 0){
    V <- object$V[1:p, 1:p, drop = FALSE]
    colnames(V) <- labels[1:p]
    rownames(V) <- labels[1:p]
  }else{
    V <- object$V
    colnames(V) <- labels
    rownames(V) <- labels
  }
  return(V)
}


#---------------------------------------------
#' Parameters estimates of a survstan model
#'
#' @aliases estimates
#' @export
#' @param object an object of the class survstan.
#' @param ... further arguments passed to or from other methods.
#' @return  the parameters estimates of a given survstan model.
#' @examples
#' \donttest{
#' library(survstan)
#' fit <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
#' estimates(fit)
#' }
#'
estimates <- function(object, ...){
  return(object$estimates)
}


#---------------------------------------------
#' Estimated regression coefficients
#'
#' @aliases coef.survstan
#' @export
#' @param object an object of the class survstan
#' @param ... further arguments passed to or from other methods
#' @return  the estimated regression coefficients
#' @examples
#' \donttest{
#' library(survstan)
#' fit <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
#' coef(fit)
#' }
#'
coef.survstan <- function(object, ...){
  p <- object$p
  if(p>0){
    coeffs <- object$estimates[1:p]
  }else{
    coeffs <- NULL
    warning("This function supports models with at least one covariate!")
  }
  return(coeffs)
}

#---------------------------------------------
#' Confidence intervals for the regression coefficients
#'
#' @aliases confint.survstan
#' @export
#' @param object an object of the class survstan.
#' @param level the confidence level required.
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param ... further arguments passed to or from other methods.
#' @return  100(1-alpha) confidence intervals for the regression coefficients.
#' @examples
#' \donttest{
#' library(survstan)
#' fit <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
#' confint(fit)
#' }
#'
confint.survstan <- function(object, parm = NULL, level=0.95, ...){
  alpha <- 1 - level
  k <- length(object$estimates)
  p <- object$p

  conf_labels <- round(100*(c(alpha/2, 1-alpha/2)),1)
  conf_labels <- paste0(conf_labels, "%")

  parameter = names(object$estimates)
  estimate = object$estimates
  se = sqrt(diag(object$V))

  ztab <- stats::qnorm(alpha/2, lower.tail = FALSE)
  CI <- data.frame(
    lwr = estimate - ztab*se,
    upr = estimate + ztab*se
  )
  names(CI) <-  conf_labels

  if(is.null(parm)){
    return(CI)
  }else{
    include <- rownames(CI) %in% parm
    return(CI[include,])
  }
}
