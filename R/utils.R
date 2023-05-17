


reparametrization <- function(object, baseline, labels, tau, p, ...){
  estimates <- object$par
  V <- solve(-object$hessian)

  if(baseline == "exponential"){
    labels <- c(labels, "lambda")
    estimates[p+1] <- estimates[p+1]/tau
    v <- diag(p+1)
    v[p+1, p+1] <- 1/tau
    V <- v%*%V%*%v
  }else{
    if(baseline == "weibull"){
      labels <- c(labels, "alpha", "gamma")
      estimates[p+2] <- estimates[p+2]*tau
      v <- diag(p+2)
      #v[p+1, p+1] <- estimates[p+1]
      v[p+2, p+2] <- tau
      V <- v%*%V%*%v
    }else if(baseline == "lognormal"){
      labels <- c(labels, "mu", "sigma")
      estimates[p+1] <- estimates[p+1] + log(tau)
      v <- diag(p+2)
    }else{ # loglogistic
      labels <- c(labels, "alpha", "gamma")
      estimates[p+2] <- estimates[p+2]*tau
      v <- diag(p+2)
      v[p+2, p+2] <- 1/tau
      V <- v%*%V%*%v
    }
  }
  names(estimates) <- labels
  rownames(V) <- labels
  colnames(V) <- labels
  res <- list(estimates=estimates, V=V)
  return(res)
}


#---------------------------------------------
#' Extract AIC from a Fitted Model
#'
#' @aliases extractAIC.survstan
#' @description Computes the (generalized) Akaike An Information Criterion for a fitted parametric model.
#' @importFrom stats extractAIC
#' @export
#' @param fit a fitted model of the class survstan
#' @param scale optional numeric specifying the scale parameter of the model. Currently only used in the "lm" method, where scale specifies the estimate of the error variance, and scale = 0 indicates that it is to be estimated by maximum likelihood.
#' @param k numeric specifying the ‘weight’ of the equivalent degrees of freedom part in the AIC formula.
#' @param ... further arguments passed to or from other methods.
#' @return  the ANOVA table.
#' @examples
#' \donttest{
#' library(survstan)
#' fit1 <- aftreg(Surv(futime, fustat) ~ 1, data = ovarian, baseline = "weibull", init = 0)
#' fit2 <- aftreg(Surv(futime, fustat) ~ rx, data = ovarian, baseline = "weibull", init = 0)
#' fit3 <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
#' extractAIC(fit1)
#' extractAIC(fit2)
#' extractAIC(fit3)
#' }
#'

extractAIC.survstan <- function(fit, scale, k=2, ...){
  edf <- length(fit$estimates)
  loglik <- fit$loglik
  c(edf, -2*loglik + k*edf)
}
