


#---------------------------------------------
#' Model.matrix method for survstan models
#'
#' @aliases model.matrix.survstan
#' @description Reconstruct the model matrix for a survstan model.
#' @export
#' @param object an object of the class survstan.
#' @param ... further arguments passed to or from other methods.
#' @return  The model matrix (or matrices) for the fit.
#' @examples
#' \donttest{
#' library(survstan)
#' fit <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
#' model.matrix(fit)
#' }
#'
model.matrix.survstan <- function(object, ...){
  mf <- stats::model.frame(object)
  formula <- object$formula
  X <- stats::model.matrix(formula, data = mf)
  attrX <- attributes(X)
  X <- X[, -1, drop = FALSE]
  attrX$dim[2] <- ncol(X)
  attrX$assign <- attrX$assign[-1]
  attrX$dimnames[[2]] <- attrX$dimnames[[2]][-1]
  attributes(X) <- attrX
  return(X)
}
