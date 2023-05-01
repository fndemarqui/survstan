


reparametrization <- function(object, baseline, labels, tau, p, ...){
  estimates <- object$par
  V <- solve(-object$hessian)

  if(baseline == "exponential"){
    labels <- c(labels, "lambda")
    estimates[p+1] <- estimates[p+1]/tau
    v <- diag(p+1)
    v[p+1, p+1] <- 1/tau
    V <- v%*%V%*%v
  }else{
    if(baseline == "weibull"){
      labels <- c(labels, "alpha", "gamma")
      estimates[p+2] <- estimates[p+2]*tau
      v <- diag(p+2)
      #v[p+1, p+1] <- estimates[p+1]
      v[p+2, p+2] <- tau
      V <- v%*%V%*%v
    }else if(baseline == "lognormal"){
      labels <- c(labels, "mu", "sigma")
      estimates[p+1] <- estimates[p+1] + log(tau)
      v <- diag(p+2)
    }else{ # loglogistic
      labels <- c(labels, "alpha", "gamma")
      estimates[p+2] <- estimates[p+2]*tau
      v <- diag(p+2)
      v[p+2, p+2] <- 1/tau
      V <- v%*%V%*%v
    }
  }
  names(estimates) <- labels
  rownames(V) <- labels
  colnames(V) <- labels
  res <- list(estimates=estimates, V=V)
  return(res)
}




