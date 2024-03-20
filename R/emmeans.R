#' Support Functions for \pkg{emmeans}
#'
#' Functions required for compatibility of \pkg{survstan} with \pkg{emmeans}.
#' Users are not required to call these functions themselves. Instead,
#' they will be called automatically by the \code{emmeans} function
#' of the \pkg{emmeans} package.
#'
#' @name emmeans-survstan-helpers



#' @rdname emmeans-survstan-helpers
#' @param object An object of the same class as is supported by a new method.
#' @param ... Additional parameters that may be supported by the method.
recover_data.survstan <- function (object, ...){
  fcall = object$call
  frame = object$model
  emmeans::recover_data(fcall, stats::delete.response(stats::terms(object)), object$na.action,
               frame = frame, pwts = stats::weights(object), ...)
}

emm_basis.survstan = function(object, trms, xlev, grid, ...) {
  m = stats::model.frame(trms, grid, na.action = stats::na.pass, xlev = xlev)
  X = stats::model.matrix(trms, m, contrasts.arg = object$contrasts)[, -1, drop = FALSE]
  bhat = stats::coef(object)
  Xmat = stats::model.matrix(trms, data=object$mf)[, -1, drop = FALSE]                      # 5
  V <- vcov(object)
  nbasis = matrix(NA)
  dfargs = list(df = Inf)
  dffun = function(k, dfargs) dfargs$df
  list(X = X, bhat = bhat, nbasis = nbasis, V = V,                  #10
       dffun = dffun, dfargs = dfargs)
}


#' @rdname emmeans-survstan-helpers
#' @param object An object of the same class as is supported by a new method.
#' @param term character specifying whether short or long term regression coefficients are to be used.
#' @param ... Additional parameters that may be supported by the method.
recover_data.ypreg <- function (object, term = c("short", "long"), ...){
  term <- match.arg(term)
  frame <- object$model
  fcall = object$call
  if (term %in% c("short", "long"))
    trms = stats::delete.response(stats::terms(object, model = term))
  emmeans::recover_data(fcall, trms, object$na.action,
               frame = frame, pwts = stats::weights(object), ...)
}

emm_basis.ypreg = function(object, trms, xlev, grid, term = c("short", "long"), ...){
  term <- match.arg(term)
  m = stats::model.frame(trms, grid, na.action = stats::na.pass, xlev = xlev)
  X = stats::model.matrix(trms, m, contrasts.arg = object$contrasts)[, -1, drop = FALSE]
  bhat = stats::coef(object)
  p <- length(bhat)/2
  if(term=="short"){
    bhat <- bhat[1:p]
    V <- vcov(object)[1:p, 1:p]
  }else{
    bhat <- bhat[(p+1):(2*p)]
    V <- vcov(object)[(p+1):(2*p), (p+1):(2*p)]
  }

  Xmat = stats::model.matrix(trms, data=object$mf)[, -1, drop = FALSE]
  nbasis = matrix(NA)
  dfargs = list(df = Inf)
  dffun = function(k, dfargs) dfargs$df
  list(X = X, bhat = bhat, nbasis = nbasis, V = V,
       dffun = dffun, dfargs = dfargs)
}


#' @rdname emmeans-survstan-helpers
#' @param object An object of the same class as is supported by a new method.
#' @param term character specifying whether AF or RH term regression coefficients are to be used.
#' @param ... Additional parameters that may be supported by the method.
recover_data.ehreg <- function (object, term = c("AF", "RH"), ...){
  term <- match.arg(term)
  frame <- object$model
  fcall = object$call
  if (term %in% c("AF", "RH"))
    trms = stats::delete.response(stats::terms(object, model = term))
  emmeans::recover_data(fcall, trms, object$na.action,
                        frame = frame, pwts = stats::weights(object), ...)
}

emm_basis.ehreg = function(object, trms, xlev, grid, term = c("AF", "RH"), ...){
  term <- match.arg(term)
  m = stats::model.frame(trms, grid, na.action = stats::na.pass, xlev = xlev)
  X = stats::model.matrix(trms, m, contrasts.arg = object$contrasts)[, -1, drop = FALSE]
  bhat = stats::coef(object)
  p <- length(bhat)/2
  if(term=="AF"){
    bhat <- bhat[1:p]
    V <- vcov(object)[1:p, 1:p]
  }else{
    bhat <- bhat[(p+1):(2*p)]
    V <- vcov(object)[(p+1):(2*p), (p+1):(2*p)]
  }

  Xmat = stats::model.matrix(trms, data=object$mf)[, -1, drop = FALSE]
  nbasis = matrix(NA)
  dfargs = list(df = Inf)
  dffun = function(k, dfargs) dfargs$df
  list(X = X, bhat = bhat, nbasis = nbasis, V = V,
       dffun = dffun, dfargs = dfargs)
}


