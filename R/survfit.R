
surv_aft <- function(time, pars, lp, baseline, p){
  time <- time*exp(-lp)
  surv <- switch(baseline,
               exponential = stats::pexp(time, rate = pars[p+1], lower.tail = FALSE),
               weibull = stats::pweibull(time, shape = pars[p+1], scale = pars[p+2], lower.tail = FALSE),
               lognormal = stats::plnorm(time, meanlog = pars[p+1], sdlog = pars[p+2], lower.tail = FALSE),
               loglogistic = actuar::pllogis(time, shape = pars[p+1], scale = pars[p+2], lower.tail = FALSE)
  )
  return(surv)
}

surv_ah <- function(time, pars, lp, baseline, p){
  time <- time*exp(lp)
  H0 <- switch(baseline,
                 exponential = -stats::pexp(time, rate = pars[p+1], lower.tail = FALSE, log.p = TRUE),
                 weibull = -stats::pweibull(time, shape = pars[p+1], scale = pars[p+2], lower.tail = FALSE, log.p = TRUE),
                 lognormal = -stats::plnorm(time, meanlog = pars[p+1], sdlog = pars[p+2], lower.tail = FALSE, log.p = TRUE),
                 loglogistic = -actuar::pllogis(time, shape = pars[p+1], scale = pars[p+2], lower.tail = FALSE, log.p = TRUE)
  )
  surv <- exp(-H0*exp(-lp))
  return(surv)
}

surv_ph <- function(time, pars, lp, baseline, p){
  H0 <- switch(baseline,
               exponential = -stats::pexp(time, rate = pars[p+1], lower.tail = FALSE, log.p = TRUE),
               weibull = -stats::pweibull(time, shape = pars[p+1], scale = pars[p+2], lower.tail = FALSE, log.p = TRUE),
               lognormal = -stats::plnorm(time, meanlog = pars[p+1], sdlog = pars[p+2], lower.tail = FALSE, log.p = TRUE),
               loglogistic = -actuar::pllogis(time, shape = pars[p+1], scale = pars[p+2], lower.tail = FALSE, log.p = TRUE)
  )
  surv <- exp(-H0*exp(lp))
  return(surv)
}

surv_po <- function(time, pars, lp, baseline, p){
  H0 <- switch(baseline,
               exponential = -stats::pexp(time, rate = pars[p+1], lower.tail = FALSE, log.p = TRUE),
               weibull = -stats::pweibull(time, shape = pars[p+1], scale = pars[p+2], lower.tail = FALSE, log.p = TRUE),
               lognormal = -stats::plnorm(time, meanlog = pars[p+1], sdlog = pars[p+2], lower.tail = FALSE, log.p = TRUE),
               loglogistic = -actuar::pllogis(time, shape = pars[p+1], scale = pars[p+2], lower.tail = FALSE, log.p = TRUE)
  )
  Rt = expm1(H0)*exp(lp)
  surv <- exp(-log1p(Rt))
  return(surv)
}


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


#---------------------------------------------
#' survfit method for survstan models
#'
#' @aliases survfit.survstan
#' @description Computes the predicted survivor function for a phpe model.
#' @importFrom survival survfit
#' @export
#' @param formula an object of the class survstan
#' @param newdata a data frame containing the set of explanatory variables.
#' @param ... further arguments passed to or from other methods.
#' @return  a list containing the estimated survival probabilities.
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
survfit.survstan <- function(formula, newdata, ...){
  object <- formula
  baseline <- object$baseline
  survreg <- object$survreg
  mf <- object$mf
  Terms <- stats::delete.response(object$terms)
  labels <- names(mf)[-1]
  time <- c(0,sort( stats::model.response(mf)[,1]))
  labels <- match.arg(names(newdata), labels, several.ok = TRUE)
  X <- stats::model.matrix(Terms, data = newdata)[, -1, drop = FALSE]
  pars <- estimates(object)
  p <- object$p
  beta <- pars[1:p]
  lp <- as.numeric(X%*%beta)
  if(survreg == "yp"){
    phi <- pars[(p+1):(2*p)]
    lp_short <- lp
    lp_long <- as.numeric(X%*%phi)
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
    "aft" = with(df, surv_aft(time, pars, lp, baseline, p)),
    "ph" = with(df, surv_ph(time, pars, lp, baseline, p)),
    "po" = with(df, surv_po(time, pars, lp, baseline, p)),
    "ah" = with(df, surv_ah(time, pars, lp, baseline, p)),
    "yp" = with(df, surv_yp(time, pars, lp_short, lp_long, baseline, p))
  )

  surv <- df %>%
    dplyr::mutate(
      id = as.factor(rep(1:J, n)),
      surv  = surv
    ) %>%
    dplyr::relocate(.data$id, .before = time)

  return(surv)
}




