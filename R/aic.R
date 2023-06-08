# internal function used to compute AIC (several models)
get_arg_names <- function(...) {
  argnames <- sys.call()
  paste0(lapply(argnames[-1], as.character))
}

#' Akaike information criterion
#' @aliases AIC.survstan
#' @export
#' @param object an object of the class survstan.
#' @param ... further arguments passed to or from other methods.
#' @param k numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.
#' @return  the Akaike information criterion value when a single model is passed to the function; otherwise, a data.frame with the Akaike information criterion values and the number of parameters is returned.
#' @examples
#' \donttest{
#' library(survstan)
#' fit1 <- aftreg(Surv(futime, fustat) ~ 1, data = ovarian, baseline = "weibull", init = 0)
#' fit2 <- aftreg(Surv(futime, fustat) ~ rx, data = ovarian, baseline = "weibull", init = 0)
#' fit3 <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
#' AIC(fit1, fit2, fit3)
#' }
#'
AIC.survstan <- function(object, ..., k = 2){
  objects <- c(as.list(environment()), list(...))
  argnames <- sys.call()
  argnames <- paste0(lapply(argnames[-1], as.character))
  k <- objects[[2]]
  objects <- objects[-2]
  J <- nargs()
  aic <- c()
  npars <- c()
  for(j in 1:J){
    loglik <- objects[[j]]$loglik
    npars[j] <- length(objects[[j]]$estimates)
    aic[j] <- -2*loglik + k*npars[j]
  }
  if(length(argnames)>1){
    # names(aic) <- argnames
    # aic <- sort(aic)
    aic <- data.frame(
      fit = argnames,
      aic = aic,
      npars = npars
    ) %>%
      dplyr::arrange(aic)
  }
  return(aic)
}
