

surv_yp <- function(time, pars, lp_short, lp_long, baseline, p){
  p <- 2*p
  H0 <- switch(baseline,
               exponential = -stats::pexp(time, rate = pars[p+1], lower.tail = FALSE, log.p = TRUE),
               weibull = -stats::pweibull(time, shape = pars[p+1], scale = pars[p+2], lower.tail = FALSE, log.p = TRUE),
               lognormal = -stats::plnorm(time, meanlog = pars[p+1], sdlog = pars[p+2], lower.tail = FALSE, log.p = TRUE),
               loglogistic = -actuar::pllogis(time, shape = pars[p+1], scale = pars[p+2], lower.tail = FALSE, log.p = TRUE)
  )
  ratio <- exp(lp_short - lp_long)
  Rt = expm1(H0)*ratio
  theta <- exp(lp_long)
  surv <- exp(-theta*log1p(Rt))
  return(surv)
}


diffSurv <- function(time, X1, X2, pars, baseline, p){
    lp_short1 <- as.numeric(X1%*%pars[1:p])
    lp_short2 <- as.numeric(X2%*%pars[1:p])
    lp_long1 <- as.numeric(X1%*%pars[(p+1):(2*p)])
    lp_long2 <- as.numeric(X2%*%pars[(p+1):(2*p)])
    St1 <- surv_yp(time, pars, lp_short1, lp_long1, baseline, p)
    St2 <- surv_yp(time, pars, lp_short2, lp_long2, baseline, p)
  return(St1-St2)
}

ypregCrossSurv <- function(X1, X2, tau0=tau0, tau=tau, pars, baseline, p){
  I <- c(tau0, 1.5*tau)
  t <- try(stats::uniroot(diffSurv, interval=I, X1=X1, X2=X2, pars=pars, baseline=baseline, p=p)$root, TRUE)
  if(is(t, "try-error")){
    return(NA)
  }else{
    return(t)
  }
}



#---------------------------------------------
#' Generic S3 method cross_time
#' @aliases cross_time
#' @export
#' @param object a fitted model object
#' @param ... further arguments passed to or from other methods.
#' @return the crossing survival time
#'
cross_time <- function(object, ...) UseMethod("cross_time")


#' Computes the crossing survival times
#'
#' @aliases cross_time.ypreg
#' @rdname cross_time-methods
#' @method cross_time ypreg
#' @export
#' @export cross_time
#' @param object an object of class ypreg
#' @param newdata1 a data frame containing the first set of explanatory variables
#' @param newdata2 a data frame containing the second set of explanatory variables
#' @param conf.level level of the confidence/credible intervals
#' @param nboot number of bootstrap samples (default nboot=1000).
#' @param ... further arguments passed to or from other methods.
#' @return  the crossing survival time
#' @examples
#' \donttest{
#' library(survstan)
#' data(ipass)
#' fit <- ypreg(Surv(time, status)~arm, data=ipass, baseline = "weibull")
#' summary(fit)
#' newdata1 <- data.frame(arm=0)
#' newdata2 <- data.frame(arm=1)
#' tcross <- cross_time(fit, newdata1, newdata2, nboot = 10)
#' tcross
#' }
#'
cross_time.ypreg <- function(object, newdata1, newdata2,
                           conf.level=0.95, nboot=1000, ...){
  baseline <- object$baseline
  survreg <- object$survreg
  mf <- object$mf
  Terms <- stats::delete.response(object$terms)
  labels <- names(mf)[-1]
  time <- sort( stats::model.response(mf)[,1])
  labels <- match.arg(names(newdata1), names(newdata2), several.ok=TRUE)
  labels <- match.arg(names(mf)[-1], names(newdata1), several.ok=TRUE)
  X1 <- stats::model.matrix(Terms, data = newdata1)[, -1, drop = FALSE]
  X2 <- stats::model.matrix(Terms, data = newdata2)[, -1, drop = FALSE]
  pars <- estimates(object)
  p <- object$p
  beta <- pars[1:p]
  phi <- pars[(p+1):(2*p)]

  lp_short1 <- as.numeric(X1%*%beta)
  lp_short2 <- as.numeric(X2%*%beta)
  lp_long1 <- as.numeric(X1%*%phi)
  lp_long2 <- as.numeric(X2%*%phi)

  tau0 <- min(time)
  tau <- object$tau

  alpha <- 1 - conf.level
  prob <- c(alpha/2, 1-alpha/2)

  t <- c()

  for(i in 1:nrow(newdata1)){
    t[i] <- ypregCrossSurv(X1=X1[i,], X2=X2[i,], tau0=tau0, tau=tau, pars=pars, baseline, p)
  }
  pars <- ypreg_boot(object, nboot=nboot)

  ci <- matrix(nrow=nrow(newdata1), ncol=2)
  for(i in 1:nrow(newdata1)){
    aux <- apply(pars, 1, ypregCrossSurv, X1=X1[i,], X2=X2[i,], tau0=tau0, tau=tau, baseline=baseline, p=p)
    ci[i,] <- stats::quantile(aux, probs=prob, na.rm=TRUE)
  }

  t <- data.frame(cbind(t, ci))
  names(t) <- c("Est.", paste(100*prob, "%", sep=""))
  return(t)
}

