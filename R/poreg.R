
#' Fitting Proportional Odds Models
#' @aliases poreg
#' @export
#' @description Function to fit proportional odds (PO) models.
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which function is called.
#' @param baseline the chosen baseline distribution; options currently available are: exponential, weibull, lognormal, loglogistic and Birnbaum-Saunders (fatigue) distributions.
#' @param dist alternative way to specify the baseline distribution (for compatibility with the \code{\link[survival]{survreg}} function); default is NULL.
#' @param init initial values specification (default value is 0); see the detailed documentation for \code{init} in \code{\link[rstan]{optimizing}}.
#' @param ... further arguments passed to other methods.
#' @return poreg returns an object of class "poreg" containing the fitted model.
#' @examples
#' \donttest{
#' library(survstan)
#' fit <- poreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull")
#' summary(fit)
#' }
#'
#'
poreg <- function(formula, data, baseline = "weibull", dist = NULL, init = 0, ...){
  if(!is.null(dist)){
    baseline <- dist
  }
  baseline <- tolower(baseline)
  baseline <- match.arg(baseline, survstan_distributions)
  Call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  resp <- stats::model.response(mf)
  X <- stats::model.matrix(mt, mf)
  time <- resp[,1]
  event <- resp[, 2]
  labels <- colnames(X)[-1]
  X <- matrix(X[,-1, drop = FALSE], ncol=length(labels))
  n <- length(time)
  p <- ncol(X)
  tau <- max(time)
  y <- time/tau
  offset <- stats::model.offset(mf)
  if(is.null(offset)){
    offset <- rep(0, n)
  }

  output <- list(call = Call, formula = stats::formula(mt), offset = offset,
                 terms = mt, mf = mf, baseline = baseline, survreg = "po",
                 n = n, p = p, tau = tau, labels = labels)

  if(init == 0 & baseline == "ggprentice"){
    init <- inits("po", p)
  }
  baseline <- set_baseline(baseline)

  stan_data <- list(time=y, event=event, X=X, n=n, p=p, offset = offset,
                    baseline=baseline, survreg = 3, tau = tau)
  fit <- rstan::optimizing(stanmodels$survreg, data = stan_data, hessian = TRUE, init = init, ...)
  res <- reparametrization(fit, survreg = "po", output$baseline, labels, tau, p)
  output$estimates <- res$estimates
  output$V <- res$V
  output$loglik = fit$value
  output$return_code = fit$return_code

  pars <- output$estimates
  if(p==0){
    lp <- 0 + offset
  }else{
    lp <- as.numeric(X%*%pars[1:p]) + offset
  }

  H0 <- cumhaz(time, pars, baseline, p)
  Rt <- expm1(H0)*exp(lp)
  output$residuals <- log1p(Rt)
  output$event <- event


  class(output) <- c("poreg", "survstan")
  return(output)
}
