
#' Fitting Proportional Hazards Models
#' @aliases phreg
#' @export
#' @description Function to fit proportional hazards (PH) models.
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which function is called.
#' @param baseline the chosen baseline distribution; options currently available are: exponential, weibull, lognormal and loglogistic distributions.
#' @param dist alternative way to specify the baseline distribution (for compability with the \code{\link[survival]{survreg}} function); default is NULL.
#' @param init initial values specification (default value is 0); see the detailed documentation for \code{init} in \code{\link[rstan]{optimizing}}.
#' @param ... further arguments passed to other methods.
#' @return phreg returns an object of class "phreg" containing the fitted model.
#' @examples
#' \donttest{
#' library(survstan)
#' fit <- phreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
#' summary(fit)
#' }
#'
#'
phreg <- function(formula, data, baseline = c("exponential", "weibull", "lognormal", "loglogistic"), dist = NULL, init = 0, ...){
  if(!is.null(dist)){
    baseline <- dist
  }
  baseline <- tolower(baseline)
  baseline <- match.arg(baseline)
  mf <- stats::model.frame(formula, data)
  resp <- stats::model.response(mf)
  time <- resp[,1]
  event <- resp[, 2]
  X <- stats::model.matrix(formula, data = mf, rhs = 1)
  labels <- colnames(X)[-1]
  X <- matrix(X[,-1, drop = FALSE], ncol=length(labels))
  n <- length(time)
  p <- ncol(X)
  tau <- max(time)
  y <- time/tau

  output <- list(call = match.call(), formula = formula,
                 terms = stats::terms.formula(formula), mf = mf,
                 baseline=baseline, survreg = "ph",
                 n=n, p=p, tau=tau, labels=labels)

  baseline <- switch(baseline,
    "exponential" = 1,
    "weibull" = 2,
    "lognormal" = 3,
    "loglogistic" = 4
  )

  stan_data <- list(time=y, event=event, X=X, n=n, p=p,
                    baseline=baseline, survreg = 2, tau = tau)
  fit <- rstan::optimizing(stanmodels$survreg, data = stan_data, hessian = TRUE, init = init, ...)
  res <- reparametrization(fit, survreg = "ph", output$baseline, labels, tau, p)
  output$estimates <- res$estimates
  output$V <- res$V
  output$loglik = fit$value
  output$return_code = fit$return_code

  pars <- output$estimates
  if(p==0){
    lp <- 0
  }else{
    lp <- as.numeric(X%*%pars[1:p])
  }

  H0 <- switch(output$baseline,
               exponential = -stats::pexp(time, rate = pars[p+1], lower.tail = FALSE, log.p = TRUE),
               weibull = -stats::pweibull(time, shape = pars[p+1], scale = pars[p+2], lower.tail = FALSE, log.p = TRUE),
               lognormal = -stats::plnorm(time, meanlog = pars[p+1], sdlog = pars[p+2], lower.tail = FALSE, log.p = TRUE),
               loglogistic = -actuar::pllogis(time, shape = pars[p+1], scale = pars[p+2], lower.tail = FALSE, log.p = TRUE)
  )
  output$residuals <- H0*exp(lp)
  output$event <- event

  class(output) <- c("phreg", "survstan")
  return(output)
}
