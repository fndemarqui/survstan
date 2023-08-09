


reparametrization <- function(object, survreg, baseline, labels, tau, p, ...){
  estimates <- object$par
  V <- MASS::ginv(-object$hessian)

  if(survreg == "yp"){
    p <- 2*p
    labels <- c(paste0("short-", labels), paste0("long-", labels))
  }

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
    }else if(baseline == "loglogistic"){ # loglogistic
      labels <- c(labels, "alpha", "gamma")
      estimates[p+2] <- estimates[p+2]*tau
      v <- diag(p+2)
      v[p+2, p+2] <- 1/tau
      V <- v%*%V%*%v
    }else if(baseline == "fatigue"){ # fatigue
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


#---------------------------------------------
#' Extract Log-Likelihood from a Fitted Model
#'
#' @aliases logLik.survstan
#' @description Extracts the log-likelihood function for a fitted parametric model.
#' @importFrom stats logLik
#' @export
#' @param object a fitted model of the class survstan
#' @param ... further arguments passed to or from other methods.
#' @return  the log-likelihood value when a single model is passed to the function; otherwise, a data.frame with the log-likelihood values and the number of parameters is returned.
#' @examples
#' \donttest{
#' library(survstan)
#' fit1 <- aftreg(Surv(futime, fustat) ~ 1, data = ovarian, baseline = "weibull", init = 0)
#' fit2 <- aftreg(Surv(futime, fustat) ~ rx, data = ovarian, baseline = "weibull", init = 0)
#' fit3 <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
#' logLik(fit1, fit2, fit3)
#' }
#'

logLik.survstan <- function(object, ...){
  objects <- c(as.list(environment()), list(...))
  argnames <- sys.call()
  argnames <- paste0(lapply(argnames[-1], as.character))
  J <- nargs()
  loglik <- c()
  npars <- c()
  for(j in 1:J){
    loglik[j] <- objects[[j]]$loglik
    npars[j] <- length(objects[[j]]$estimates)
  }
  if(length(argnames)>1){
    # names(loglik) <- argnames
    # loglik <- sort(loglik, decreasing = TRUE)
    loglik <- data.frame(
      fit = argnames,
      loglik = loglik,
      npars = npars
    ) %>%
      dplyr::arrange(desc(loglik))
  }
  return(loglik)
}


survstan_distributions <- c("exponential", "weibull", "lognormal", "loglogistic", "fatigue")
