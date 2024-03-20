#' Print the summary.survstan output
#'
#' @export
#' @description Produces a printed summary of a fitted survstan model.
#' @param x an object of the class summary.survstan.
#' @param ... further arguments passed to or from other methods.
#' @return No return value, called for side effects.
#'
print.summary.survstan <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat(x$message)
  cat("\n")
  if(x$p>0){
    cat("Regression coefficients:\n")
    stats::printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
    cat("\n")
    cat("Baseline parameters:\n")
  }
  stats::printCoefmat(x$tbl, P.value=FALSE, has.Pvalue=FALSE)
  cat("--- \n")
  cat("loglik =", x$loglik, " ", "AIC =", x$AIC,"\n")
}


#---------------------------------------------

#' Summary for a survstan object
#'
#' @aliases summary.survstan
#' @export
#' @param object the result of a call to summary.survstan
#' @param conf.level the confidence level required.
#' @param ... further arguments passed to or from other methods.
#' @return an object of the class summary.survstan containing a summary of the fitted model.
#'
summary.survstan <- function(object, conf.level = 0.95, ...){

  labels <- object$labels
  baseline <- object$baseline
  if(object$survreg == "yp"){
    p <- 2*object$p
    labels <- c(paste0("short-", labels), paste0("long-", labels))
  }else if(object$survreg == "eh"){
    p <- 2*object$p
    labels <- c(paste0("AF-", labels), paste0("RH-", labels))
  }else{
    p <- object$p
  }
  n <- object$n
  survreg <- object$survreg

  estimates <- object$estimates
  SE <- se(object)
  V <- object$V
  k <- length(estimates)
  tau <- object$tau

  loglik <- object$loglik
  AIC <- -2*loglik + 2*k
  BIC <- -2*loglik + k*log(n)



  if(p>0){
    coefficients <- estimates[1:p]
    SE <- SE[1:p]
    zval <- coefficients / SE
    TAB <- cbind(Estimate = coefficients,
                 `Std. Error` = SE,
                 `z value` = zval,
                 `Pr(>|z|)` = 2*stats::pnorm(-abs(zval)))

    rownames(TAB) <- labels
    estimates <- estimates[-(1:p)]
    SE <- se(object, all = TRUE)[-(1:p)]
  }


  CI <- stats::confint(object, parm = names(estimates), level = conf.level)
  lwr <- CI[, 1]
  upr <- CI[, 2]
  tbl <- cbind(estimates, SE, lwr, upr)

  alpha <- 1-conf.level
  conf_labels <- round(100*(c(alpha/2, 1-alpha/2)),1)
  conf_labels <- paste0(conf_labels, "%")
  colnames(tbl) <- c("Estimate", "Std. Error", conf_labels)

  if(baseline == "ggprentice"){
    baseline <- "generalized gamma (Prentice)"
  }else if(baseline == "ggstacy"){
    baseline <- "generalized gamma (Prentice)"
  }


  message <- switch(survreg,
    "aft" = paste0("Accelerated failure time model fit with ", baseline , " baseline distribution: \n"),
    "ph" = paste0("Proportional hazards model fit with ", baseline , " baseline distribution: \n"),
    "po" = paste0("Proportional odds model fit with ", baseline , " baseline distribution: \n"),
    "ah" = paste0("Accelerated hazard model fit with ", baseline , " baseline distribution: \n"),
    "yp" = paste0("Yang & Prentice model fit with ", baseline , " baseline distribution: \n"),
    "eh" = paste0("Extended hazard model fit with ", baseline , " baseline distribution: \n")
  )

  if(p>0){
    res <- list(call=object$call,
                coefficients=TAB,
                tbl = tbl,
                loglik=loglik, AIC=AIC, message=message, p=p)
  }else{
    message <- paste("Fit of", baseline, "distribution:")
    res <- list(call=object$call,
                tbl = tbl,
                loglik=loglik, AIC=AIC, message=message, p=p)

  }

  class(res) <- "summary.survstan"

  return(res)
}


