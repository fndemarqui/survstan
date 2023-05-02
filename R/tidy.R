
#---------------------------------------------
#' Generic S3 method tidy
#' @aliases tidy
#' @export
#' @param object a fitted model object.
#' @param conf.level the confidence level required.
#' @details Convert a fitted model into a tibble.
#' @param ... further arguments passed to or from other methods.
#' @return a tibble with a summary of the fit.
#'

tidy <- function(object, conf.level = 0.95, ...) UseMethod("tidy")

#' Tidy a survstan object
#' @aliases tidy.survstan
#' @export
#' @param object a fitted model object.
#' @param conf.level the confidence level required.
#' @details Convert a fitted model into a tibble.
#' @param ... further arguments passed to or from other methods.
#' @return a tibble with a summary of the fit.
#' @examples
#' \donttest{
#' library(survstan)
#' fit <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
#' tidy(fit)
#' }
#'
tidy.survstan <- function(object, conf.level = 0.95, ...){
  alpha <- 1 - conf.level
  k <- length(object$estimates)
  p <- object$p

  conf_labels <- round(100*(c(alpha/2, 1-alpha/2)),1)
  conf_labels <- paste0(conf_labels, "%")

  parameter = names(object$estimates)
  estimate = object$estimates
  se = sqrt(diag(object$V))

  ztab <- stats::qnorm(alpha/2, lower.tail = FALSE)
  tbl <- tibble::tibble(
    parameter = parameter,
    type = c(rep("coefficient", p), rep("baseline", k-p)),
    estimate = estimate,
    se = se,
    lwr = estimate - ztab*se,
    upr = estimate + ztab*se,
  )

  names(tbl) <- c(names(tbl)[1:4], conf_labels)
  return(tbl)

}
