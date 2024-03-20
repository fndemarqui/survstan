#' The Gompertz Distribution
#' @name Gompertz
#' @aliases gompertz
#' @aliases Gompertz dgompertz pgompertz qgompertz hgompertz Hgompertz
#'
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#' @description Probability function, distribution function, quantile function and random generation for the  distribution with parameters alpha and gamma.
#'
#' @param x	vector of (non-negative integer) quantiles.
#' @param q	vector of quantiles.
#' @param p	vector of probabilities.
#' @param n	number of random values to return.
#' @param alpha shape parameter of the  distribution (alpha > 0).
#' @param gamma scale parameter of the  distribution (gamma > 0).
#' @param log,log.p	 logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	logical; if TRUE (default), probabilities are \eqn{P[X \le x]}; otherwise, \eqn{P[X > x]}.
#' @param ... further arguments passed to other methods.
#'
#' @details
#'
#' Probability density function:
#' \deqn{
#' f(x|\alpha, \gamma) = \alpha\gamma \exp\{\gamma x - \alpha(e^{\gamma x} - 1)\}I_{[0, \infty)}(x),
#' }
#' for \eqn{\alpha>0} and \eqn{\gamma>0}.
#'
#' Distribution function:
#' \deqn{
#' F(x|\alpha, \gamma) = 1 - \exp\{- \alpha(e^{\gamma x} - 1)\},
#' }
#' for \eqn{x>0}, \eqn{\alpha>0} and \eqn{\gamma>0}.
#' @return dgompertz gives the (log) probability function, pgompertz gives the (log) distribution function, qgompertz gives the quantile function, and rgompertz generates random deviates.
#'

##' @export
##' @rdname Gompertz
dgompertz <- function(x, alpha = 1, gamma = 1,  log = FALSE, ...){
  lf <- log(alpha) + log(gamma) + gamma*x - alpha*expm1(gamma*x)
  if(isTRUE(log)){
    return(lf)
  }else{
    return(exp(lf))
  }
}

##' @export
##' @rdname Gompertz
pgompertz <- function(q, alpha = 1, gamma = 1,  lower.tail = TRUE, log.p = FALSE, ...){
  aux <- -alpha*expm1(gamma*q)
  if(isTRUE(lower.tail)){
    p <- -expm1(aux)
  }else{
    p <- exp(aux)
  }
  if(isTRUE(log.p)){
    return(log(p))
  }else{
    return(p)
  }
}

##' @export
##' @rdname Gompertz
qgompertz <- function(p, alpha = 1, gamma = 1, lower.tail = FALSE, log.p = FALSE, ...){
  if(isTRUE(lower.tail)){
    v <- -log(p)/alpha + 1
  }else{
    v <- -log(1-p)/alpha + 1
  }
  t <- log(v)/gamma
  return(t)
}

##' @export
##' @rdname Gompertz
rgompertz <- function(n, alpha = 1, gamma = 1, ...){
  u <- stats::runif(n)
  x <- qgompertz(u, alpha, gamma)
  return(x)
}
