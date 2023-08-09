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

  if(x$p>0){
    cat(x$message)
    cat("\n")
    cat("Regression coefficients:\n")
    stats::printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
    cat("\n")
  }

  cat("Baseline parameters:\n")
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
  }else{
    p <- object$p
  }
  n <- object$n
  survreg <- object$survreg

  estimates <- object$estimates
  V <- object$V
  k <- length(estimates)
  tau <- object$tau

  loglik <- object$loglik
  AIC <- -2*loglik + 2*k
  BIC <- -2*loglik + k*log(n)



  if(p>0){
    coefficients <- estimates[1:p]
    vcov <- V[1:p, 1:p, drop = FALSE]

    se <- sqrt(diag(vcov))
    zval <- coefficients / se
    TAB <- cbind(Estimate = coefficients,
                 StdErr = se,
                 z.value = zval,
                 p.value = 2*stats::pnorm(-abs(zval)))

    rownames(TAB) <- labels

    estimate <- estimates[-(1:p)]
    vcov <- V[-(1:p), -(1:p), drop = FALSE]
  }else{
    estimate <- estimates
    vcov <- V
  }


  # baseline parameters:
  se <- sqrt(diag(vcov))
  alpha <- 1-conf.level
  ztab <- stats::qnorm(alpha/2, lower.tail=FALSE)
  d <- ztab*se
  lwr <- estimate - d
  upr <- estimate + d
  tbl <- cbind(estimate, se, lwr, upr)

  conf_labels <- round(100*(c(alpha/2, 1-alpha/2)),1)
  conf_labels <- paste0(conf_labels, "%")
  colnames(tbl) <- c(colnames(tbl)[1:2], conf_labels)

  message <- switch(survreg,
    "aft" = paste0("Accelerated failure time model fit with ", baseline , " baseline distribution \n"),
    "ph" = paste0("Proportional hazards model fit with ", baseline , " baseline distribution \n"),
    "po" = paste0("Proportional odds model fit with ", baseline , " baseline distribution \n"),
    "ah" = paste0("Accelerated hazard model fit with ", baseline , " baseline distribution \n"),
    "yp" = paste0("Yang & Prentice model fit with ", baseline , " baseline distribution \n")
  )

  if(p>0){
    res <- list(call=object$call,
                coefficients=TAB,
                tbl = tbl,
                loglik=loglik, AIC=AIC, message=message, p=p)
  }else{
    res <- list(call=object$call,
                tbl = tbl,
                loglik=loglik, AIC=AIC, message=message, p=p)

  }

  class(res) <- "summary.survstan"

  return(res)
}


