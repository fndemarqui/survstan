

ypreg_boot <- function(object, nboot, ...){
  baseline <- object$baseline
  formula <- object$formula
  formula <- stats::update(formula, Surv(time, status) ~ .)
  mf <- object$mf
  resp <- stats::model.response(mf)
  time <- resp[,1]
  status <- resp[,2]
  data <- data.frame(cbind(time, status, mf[,-1]))
  names(data) <- c("time", "status", names(mf)[-1])
  n <- object$n
  p <- object$p
  tau <- object$tau

  index <- 1:n
  index1 <- which(status==1)
  index2 <- which(status==0)
  n1 <- length(index1)
  n2 <- length(index2)
  pars_hat <- estimates(object)
  k <- length(estimates(object))
  pars <- matrix(nrow=nboot, ncol=k)

  for(b in 1:nboot){
    samp1 <- sample(index1, size=n1, replace=TRUE)
    samp2 <- sample(index2, size=n2, replace=TRUE)
    samp <- c(samp1, samp2)
    suppressWarnings({invisible(utils::capture.output(fit <- ypreg(formula, data=data[samp,], baseline = baseline)))})
    if(!is(object, "try-error")){
      pars[b, ] <- estimates(fit)
    }

  }

  colnames(pars) <- names(pars_hat)

  return(pars)
}
