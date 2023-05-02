
#' Fitting Accelerated Failure Time Models
#' @aliases aftreg
#' @export
#' @description Function to fit accelerated failure time (AFT)  models.
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which function is called.
#' @param baseline the chosen baseline distribution; options currently available are: exponential, weibull, lognormal and loglogistic distributions.
#' @param ... further arguments passed to other methods.
#' @return aftreg returns an object of class "aftreg" containing the fitted model.
#' @examples
#' \donttest{
#' library(survstan)
#' fit <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
#' summary(fit)
#' }
#'
aftreg <- function(formula, data, baseline = c("exponential", "weibull", "lognormal", "loglogistic"), ...){
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
                 baseline=baseline, survreg = "aft",
                 n=n, p=p, tau=tau, labels=labels)

  baseline <- switch(baseline,
    "exponential" = 1,
    "weibull" = 2,
    "lognormal" = 3,
    "loglogistic" = 4
  )

  stan_data <- list(time=y, event=event, X=X, n=n, p=p,
                    baseline=baseline, survreg = 1, tau = tau)
  fit <- rstan::optimizing(stanmodels$survreg, data = stan_data, hessian = TRUE, ...)
  res <- reparametrization(fit, output$baseline, labels, tau, p)
  output$estimates <- res$estimates
  output$V <- res$V
  output$loglik = fit$value
  output$return_code = fit$return_code


  pars <- output$estimates
  if(p==0){
    nu <- time
  }else{
    lp <- as.numeric(X%*%pars[1:p])
    nu <- time*exp(-lp)
  }

  output$residuals <- switch(output$baseline,
                      exponential = -stats::pexp(nu, rate = pars[p+1], lower.tail = FALSE, log.p = TRUE),
                      weibull = -stats::pweibull(nu, shape = pars[p+1], scale = pars[p+2], lower.tail = FALSE, log.p = TRUE),
                      lognormal = -stats::plnorm(nu, meanlog = pars[p+1], sdlog = pars[p+2], lower.tail = FALSE, log.p = TRUE),
                      loglogistic = -actuar::pllogis(nu, shape = pars[p+1], scale = pars[p+2], lower.tail = FALSE, log.p = TRUE)
  )
  output$event <- event

  class(output) <- c("aftreg", "survstan")
  return(output)
}
