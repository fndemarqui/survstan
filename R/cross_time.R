
diffSurv <- function(time, X1, X2, pars, baseline, survreg, p){
    lp11 <- as.numeric(X1%*%pars[1:p])
    lp12 <- as.numeric(X2%*%pars[1:p])
    if(survreg == "yp" | survreg == "eh"){
      lp21 <- as.numeric(X1%*%pars[(p+1):(2*p)])
      lp22 <- as.numeric(X2%*%pars[(p+1):(2*p)])
    }else{
      lp11 <- 0
      lp22 <- 0
    }
    if(survreg == "yp"){
      St1 <- surv_yp(time, pars, lp11, lp21, baseline, p)
      St2 <- surv_yp(time, pars, lp12, lp22, baseline, p)
    }else if(survreg == "eh"){
      St1 <- surv_eh(time, pars, lp11, lp21, baseline, p)
      St2 <- surv_eh(time, pars, lp12, lp22, baseline, p)
    }else{
      St1 <- surv_ah(time, pars, lp11, baseline, p)
      St2 <- surv_ah(time, pars, lp12, baseline, p)
    }
  return(St1-St2)
}

crossing_time <- function(X1, X2, tau0=tau0, tau=tau, pars, baseline, survreg, p){
  I <- c(tau0, 1.5*tau)
  t <- try(stats::uniroot(diffSurv, interval=I, X1=X1, X2=X2, pars=pars, baseline=baseline, survreg=survreg, p=p)$root, TRUE)
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
#' @aliases cross_time.survstan
#' @rdname cross_time-methods
#' @method cross_time survstan
#' @export
#' @export cross_time
#' @param object an object of class survstan
#' @param newdata1 a data frame containing the first set of explanatory variables
#' @param newdata2 a data frame containing the second set of explanatory variables
#' @param conf.level level of the confidence/credible intervals
#' @param nboot number of bootstrap samples (default nboot=1000).
#' @param cores number of cores to be used in the bootstrap sampling; default is 1 core;
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
cross_time.survstan <- function(object, newdata1, newdata2,
                           conf.level=0.95, nboot=1000,
                           cores = 1, ...){

  message("Please, be patient!!!")
  if(cores == 1){
    message("Bootstrap samples draw using ", cores, " core")
  }else{
    message("Bootstrap samples draw using ", cores, " cores")
  }

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
    t[i] <- crossing_time(X1=X1[i,], X2=X2[i,], tau0=tau0, tau=tau, pars=pars, baseline=baseline, survreg=survreg, p)
  }
  pars <- bootstrap(object, nboot=nboot, cores = cores)

  ci <- matrix(nrow=nrow(newdata1), ncol=2)
  for(i in 1:nrow(newdata1)){
    aux <- apply(pars, 1, crossing_time, X1=X1[i,], X2=X2[i,], tau0=tau0, tau=tau, baseline=baseline, p=p)
    ci[i,] <- stats::quantile(aux, probs=prob, na.rm=TRUE)
  }

  t <- data.frame(cbind(t, ci))
  names(t) <- c("Est.", paste(100*prob, "%", sep=""))
  return(t)
}

