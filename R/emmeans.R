

recover_data.survstan <- function (object, frame = object$model, ...){
  fcall = object$call
  recover_data(fcall, delete.response(terms(object)), object$na.action,
               frame = frame, pwts = weights(object), ...)
}

emm_basis.survstan = function(object, trms, xlev, grid, ...) {
  m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
  X = model.matrix(trms, m, contrasts.arg = object$contrasts)[, -1, drop = FALSE]
  bhat = coef(object)
  Xmat = model.matrix(trms, data=object$mf)[, -1, drop = FALSE]                      # 5
  #V = .my.vcov(object, ...)[seq_len(k), seq_len(k), drop = FALSE]
  V <- vcov(object)
  nbasis = matrix(NA)
  dfargs = list(df = Inf)
  dffun = function(k, dfargs) dfargs$df
  list(X = X, bhat = bhat, nbasis = nbasis, V = V,                  #10
       dffun = dffun, dfargs = dfargs)
}



recover_data.ypreg <- function (object, term = c("short", "long"), ...){
  term <- match.arg(term)
  frame <- object$model
  fcall = object$call
  if (term %in% c("short", "long"))
    trms = delete.response(terms(object, model = term))
  recover_data(fcall, trms, object$na.action,
               frame = frame, pwts = weights(object), ...)
}

emm_basis.ypreg = function(object, trms, xlev, grid, term = c("short", "long"), ...){
  term <- match.arg(term)
  m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
  X = model.matrix(trms, m, contrasts.arg = object$contrasts)[, -1, drop = FALSE]
  bhat = coef(object)
  p <- length(bhat)/2
  if(term=="short"){
    bhat <- bhat[1:p]
    V <- vcov(object)[1:p, 1:p]
  }else{
    bhat <- bhat[(p+1):(2*p)]
    V <- vcov(object)[(p+1):(2*p), (p+1):(2*p)]
  }

  Xmat = model.matrix(trms, data=object$mf)[, -1, drop = FALSE]
  nbasis = matrix(NA)
  dfargs = list(df = Inf)
  dffun = function(k, dfargs) dfargs$df
  list(X = X, bhat = bhat, nbasis = nbasis, V = V,
       dffun = dffun, dfargs = dfargs)
}


