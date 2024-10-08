
surv_aft <- function(time, pars, lp, baseline, p, m){
  time <- time*exp(-lp)
  H0 <- cumhaz(time, pars, baseline, p, m)
  surv <- exp(-H0)
  return(surv)
}

surv_ah <- function(time, pars, lp, baseline, p, m){
  time <- time*exp(-lp)
  H0 <- cumhaz(time, pars, baseline, p, m)
  surv <- exp(-H0*exp(lp))
  return(surv)
}

surv_ph <- function(time, pars, lp, baseline, p, m, rho){
  H0 <- cumhaz(time, pars, baseline, p, m)
  surv <- exp(-H0*exp(lp))
  return(surv)
}

surv_po <- function(time, pars, lp, baseline, p, m, rho){
  H0 <- cumhaz(time, pars, baseline, p, m, rho)
  Rt = expm1(H0)*exp(lp)
  surv <- exp(-log1p(Rt))
  return(surv)
}


surv_yp <- function(time, pars, lp_short, lp_long, baseline, p, m, rho){
  p <- 2*p
  H0 <- cumhaz(time, pars, baseline, p, m, rho)
  ratio <- exp(lp_short - lp_long)
  Rt = expm1(H0)*ratio
  theta <- exp(lp_long)
  surv <- exp(-theta*log1p(Rt))
  return(surv)
}


surv_eh <- function(time, pars, lp1, lp2, baseline, p, m){
  p <- 2*p
  time <- time/exp(lp1)
  H0 <- cumhaz(time, pars, baseline, p, m)
  surv <- exp(- H0*exp(lp1+lp2))
  return(surv)
}

#---------------------------------------------
#' survfit method for survstan models
#'
#' @aliases survfit.survstan
#' @description Computes the predicted survivor function for a phpe model.
#' @importFrom survival survfit
#' @export
#' @param formula an object of the class survstan
#' @param newdata a data frame containing the set of explanatory variables; if NULL, a data.frame with the observed failure times and their corresponding estimated baseline survivals is returned.
#' @param ... further arguments passed to or from other methods.
#' @return  a data.frame containing the estimated survival probabilities.
#' @examples
#' \donttest{
#' library(survstan)
#' library(ggplot2)
#' data(ipass)
#' ipass$arm <- as.factor(ipass$arm)
#' fit <- ypreg(Surv(time, status)~arm, data=ipass, baseline = "weibull")
#' summary(fit)
#' newdata <- data.frame(arm=as.factor(0:1))
#' surv <- survfit(fit, newdata)
#' ggplot(surv, aes(x=time, y=surv, color = arm)) +
#'   geom_line()
#' }
#'
survfit.survstan <- function(formula, newdata = NULL, ...){
  object <- formula
  baseline <- object$baseline
  survreg <- object$survreg
  mf <- object$mf
  Terms <- stats::delete.response(object$terms)
  labels <- names(mf)[-1]
  time <- c(0,sort( stats::model.response(mf)[,1]))
  pars <- estimates(object)
  m <- object$m
  rho <- object$rho

  if(is.null(newdata)){
    H0 <- cumhaz(time, pars, baseline, p=0, m, rho)
    surv <- data.frame(
      time = time,
      surv = exp(-H0)
    )
    return(surv)
  }else{
    labels <- match.arg(names(newdata), labels, several.ok = TRUE)
    X <- stats::model.matrix(Terms, data = newdata)[, -1, drop = FALSE]
    p <- object$p
    beta <- pars[1:p]
    n <- nrow(newdata)
    mf <- stats::model.frame(Terms, data = newdata)
    offset <- stats::model.offset(mf)
    if(is.null(offset)){
      offset <- rep(0, n)
    }


    lp <- as.numeric(X%*%beta) + offset
    if(survreg == "yp"){
      phi <- pars[(p+1):(2*p)]
      lp_short <- lp
      lp_long <- as.numeric(X%*%phi) + offset
    }else if(survreg == "eh"){
      phi <- pars[(p+1):(2*p)]
      lp1 <- lp
      lp2 <- as.numeric(X%*%phi) + offset
    }

    newdata$lp <- lp
    df <- data.frame(
      time = time
    ) %>%
      dplyr::cross_join(newdata)
    N <- nrow(df)
    n <- length(time)
    J <- N/n

    surv <- switch (survreg,
                    "aft" = with(df, surv_aft(time, pars, lp, baseline, p, m)),
                    "ph" = with(df, surv_ph(time, pars, lp, baseline, p, m, rho)),
                    "po" = with(df, surv_po(time, pars, lp, baseline, p, m, rho)),
                    "ah" = with(df, surv_ah(time, pars, lp, baseline, p, m)),
                    "yp" = with(df, surv_yp(time, pars, lp_short, lp_long, baseline, p, m, rho)),
                    "eh" = with(df, surv_eh(time, pars, lp1, lp2, baseline, p, m))
    )

    surv <- df %>%
      dplyr::mutate(
        id = as.factor(rep(1:J, n)),
        surv  = surv
      ) %>%
      dplyr::relocate(.data$id, .before = time) %>%
      dplyr::select(-.data$lp)

    return(surv)
  }
}




