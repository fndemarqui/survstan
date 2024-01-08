
#' @importFrom generics tidy
#' @export
generics::tidy

#' Tidy a survstan object
#' @aliases tidy.survstan
#' @importFrom broom tidy
#' @export
#' @param x a fitted model object.
#' @param conf.int Logical indicating whether or not to include a confidence interval in the tidied output. Defaults to FALSE.
#' @param conf.level the confidence level required.
#' @details Convert a fitted model into a tibble.
#' @param ... further arguments passed to or from other methods.
#' @return a tibble with a summary of the fit.
#' @examples
#' \donttest{
#' library(survstan)
#' fit <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull")
#' tidy(fit)
#' }
#'
tidy.survstan <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {

  result <- summary(x)$coefficients %>%
    tibble::as_tibble(rownames = "term")
  colnames(result) <- c("term", "estimate", "std.error", "statistic", "p.value")
  if(conf.int){
    ci <- confint(x, level = conf.level)
    names(ci) <- c("conf.low", "conf.high")
    ci <- ci %>%
      tibble::as_tibble(rownames = "term")
    result <- dplyr::left_join(result, ci, by = "term")
  }
  return(result)
}


