

# add here new baseline distributions
set_baseline <- function(baseline){
  baseline <- switch(baseline,
                     "exponential" = 1,
                     "weibull" = 2,
                     "lognormal" = 3,
                     "loglogistic" = 4,
                     "fatigue" = 5,
                     "gamma" = 6,
                     "rayleigh" = 7
  )
  return(baseline)
}


delta_method <- function(estimates, V, pars){
  labels <- names(estimates)
  is_positive <- labels %in% pars
  estimates[!is_positive] <- 1
  D <- diag(estimates) # paremeters already transformed!!!
  V <- D%*%V%*%D
  return(V)
}

reparametrization <- function(object, survreg, baseline, labels, tau, p, ...){
  estimates <- object$par
  V <- try(chol2inv(chol(-object$hessian)), TRUE)
  if(class(V)[1] == "try-error"){
    V <- MASS::ginv(-object$hessian)
  }
  npar <- length(estimates)
  v <- diag(npar)

  if(survreg == "yp"){
    p <- 2*p
    if(p>0){
      labels <- c(paste0("short-", labels), paste0("long-", labels))
    }
  }

  if(survreg == "eh"){
    p <- 2*p
    if(p>0){
      labels <- c(paste0("AF-", labels), paste0("RH-", labels))
    }
  }

  if(baseline == "exponential"){
    labels <- c(labels, "lambda")
    names(estimates) <- labels
    colnames(v) = labels
    rownames(v) = labels
    V <- delta_method(estimates, V, "lambda")
    estimates["lambda"] <- estimates["lambda"]/tau
    diag(v)["lambda"] <- 1/tau
    V <- v%*%V%*%v
  }else if(baseline == "weibull"){
      labels <- c(labels, "alpha", "gamma")
      names(estimates) <- labels
      colnames(v) = labels
      rownames(v) = labels
      V <- delta_method(estimates, V, c("alpha", "gamma"))
      estimates["gamma"] <- estimates["gamma"]*tau
      diag(v)["gamma"] <- tau
      V <- v%*%V%*%v
  }else if(baseline == "lognormal"){
    labels <- c(labels, "mu", "sigma")
    names(estimates) <- labels
    estimates["mu"] <- estimates["mu"] + log(tau)
    V <- delta_method(estimates, V, "sigma")
    colnames(V) <- labels
    rownames(V) <- labels
  }else if(baseline == "loglogistic"){ # loglogistic
    labels <- c(labels, "alpha", "gamma")
    names(estimates) <- labels
    colnames(v) = labels
    rownames(v) = labels
    V <- delta_method(estimates, V, c("alpha", "gamma"))
    estimates["gamma"] <- estimates["gamma"]*tau
    diag(v)["gamma"] <- tau
    V <- v%*%V%*%v
  }else if(baseline == "fatigue"){ # fatigue
    labels <- c(labels, "alpha", "gamma")
    names(estimates) <- labels
    colnames(v) = labels
    rownames(v) = labels
    V <- delta_method(estimates, V, c("alpha", "gamma"))
    estimates["gamma"] <- estimates["gamma"]*tau
    diag(v)["gamma"] <- tau
    V <- v%*%V%*%v
  }else if(baseline == "gamma"){ # gamma
    labels <- c(labels, "alpha", "gamma")
    names(estimates) <- labels
    colnames(v) = labels
    rownames(v) = labels
    V <- delta_method(estimates, V, c("alpha", "gamma"))
    estimates["gamma"] <- estimates["gamma"]/tau
    diag(v)["gamma"] <- 1/tau
    V <- v%*%V%*%v
  }else if(baseline == "rayleigh"){
    labels <- c(labels, "sigma")
    names(estimates) <- labels
    colnames(v) = labels
    rownames(v) = labels
    V <- delta_method(estimates, V, "sigma")
    estimates["sigma"] <- estimates["sigma"]*tau
    diag(v)["sigma"] <- tau
    V <- v%*%V%*%v
  }
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


# This internal function returns the cumulative baseline hazard function
cumhaz <- function(time, pars, baseline, p){
  H0 <- switch(baseline,
               exponential = -stats::pexp(time, rate = pars[p+1], lower.tail = FALSE, log.p = TRUE),
               weibull = -stats::pweibull(time, shape = pars[p+1], scale = pars[p+2], lower.tail = FALSE, log.p = TRUE),
               lognormal = -stats::plnorm(time, meanlog = pars[p+1], sdlog = pars[p+2], lower.tail = FALSE, log.p = TRUE),
               loglogistic = -actuar::pllogis(time, shape = pars[p+1], scale = pars[p+2], lower.tail = FALSE, log.p = TRUE),
               fatigue = -extraDistr::pfatigue(time, alpha = pars[p+1], beta = pars[p+2], mu = 0, lower.tail = FALSE, log.p = TRUE),
               gamma = -stats::pgamma(time, shape = pars[p+1], rate = pars[p+2], lower.tail = FALSE, log.p = TRUE),
               rayleigh = -extraDistr::prayleigh(time, sigma = pars[p+1], lower.tail = FALSE, log.p = TRUE)
  )
  return(H0)
}
