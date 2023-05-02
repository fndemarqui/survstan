

extract_formulas <- function(object){
  return(object$formula)
}


#---------------------------------------------
#' anova method for survstan models
#'
#' @aliases anova.survstan
#' @description Compute analysis of variance (or deviance) tables for one or more fitted model objects.
#' @importFrom stats anova
#' @export
#' @param ... further arguments passed to or from other methods.
#' @return  the ANOVA table.
#' @examples
#' \donttest{
#' library(survstan)
#' fit1 <- aftreg(Surv(futime, fustat) ~ 1, data = ovarian, baseline = "weibull", init = 0)
#' fit2 <- aftreg(Surv(futime, fustat) ~ rx, data = ovarian, baseline = "weibull", init = 0)
#' fit3 <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
#' anova(fit1, fit2, fit3)
#' }
#'
anova.survstan <- function(...){
  models <- c(as.list(environment()), list(...))

  J <- nargs()
  labels <- paste0("Model ", 1:J, ":")
  k <- c()
  df <- c()
  k[J] <- length(models[[J]]$estimates)
  LR <- c()
  p.value <- c()
  loglik <- c()
  survreg <- c()
  loglik[J] <- models[[J]]$loglik
  survreg[J] <- models[[J]]$survreg
  for(j in 1:(J-1)){
    loglik[j] <- models[[j]]$loglik
    survreg[j] <- models[[j]]$survreg
    LR[j] <- 2*(models[[J]]$loglik - models[[j]]$loglik)
    k[j] <- length(models[[j]]$estimates)
    df[j] <- k[J]-k[j]
    p.value[j] <- stats::pchisq(LR[j], df = df[j], lower.tail = FALSE)
  }

  check <- rep(survreg[1], J)
  if(isTRUE(all.equal(check, survreg))){
    tab <- cbind("loglik" = loglik[-J], LR, df, 'Pr(>Chi)' = p.value)
    aux <- matrix(c(loglik[J], NA, NA, NA), nrow = 1)
    tab <- rbind(tab, aux)
    rownames(tab) <- labels
    formulas <- sapply(models, extract_formulas)

    cat("\n")
    for(j in 1:J){
      cat("Model", j, ": ", deparse(formulas[[j]]), "\n")
    }
    cat("--- \n")
    stats::printCoefmat(tab, P.values=TRUE, has.Pvalue = TRUE, na.print = "-")
  }else{
    warning("Only models belonging to the same regression class are allowed!")
  }


}
