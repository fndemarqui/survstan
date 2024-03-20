
#' The Generalized Gamma Distribution (Stacy's original parametrization)
#' @name ggstacy
#' @aliases ggstacy dggstacy pggstacy qggstacy rggstacy
#'
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#' @description Probability function, distribution function, quantile function and random generation for the  distribution with parameters alpha, gamma and kappa.
#'
#' @param x	vector of (non-negative integer) quantiles.
#' @param q	vector of quantiles.
#' @param p	vector of probabilities.
#' @param n	number of random values to return.
#' @param alpha shape parameter of the  distribution (alpha > 0).
#' @param gamma scale parameter of the  distribution (gamma > 0).
#' @param kappa shape parameter of the  distribution (kappa > 0).
#' @param log,log.p	 logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	logical; if TRUE (default), probabilities are \eqn{P[X \le x]}; otherwise, \eqn{P[X > x]}.
#' @param ... further arguments passed to other methods.
#'
#' @details
#'
#' Probability density function:
#' \deqn{
#' f(x|\alpha, \gamma, \kappa) = \frac{\kappa}{\gamma^{\alpha}\Gamma(\alpha/\kappa)}x^{\alpha-1}\exp\left\{-\left(\frac{x}{\gamma}\right)^{\kappa}\right\}I_{[0, \infty)}(x),
#' }
#' for \eqn{\alpha>0}, \eqn{\gamma>0} and \eqn{\kappa>0}.
#'
#' Distribution function:
#' \deqn{
#' F(t|\alpha, \gamma, \kappa) = F_{G}(x|\nu, 1),
#' }
#' where \eqn{x = \displaystyle\left(\frac{t}{\gamma}\right)^\kappa}, and \eqn{F_{G}(\cdot|\nu, 1)} corresponds to the distribution function of a gamma distribution with shape parameter \eqn{\nu = \alpha/\gamma} and scale parameter equals to 1.
#' @return dggstacy gives the (log) probability function, pggstacy gives the (log) distribution function, qggstacy gives the quantile function, and rggstacy generates random deviates.
#'


##' @export
##' @rdname ggstacy
#'
dggstacy <- function(x, alpha, gamma, kappa, log = FALSE){
  lpdf <- log(kappa) - alpha*log(gamma) - lgamma(alpha/kappa) + (alpha-1)*log(x) - (x/gamma)^kappa
  if(isTRUE(log)){
    return(lpdf)
  }else{
    return(exp(lpdf))
  }
}

##' @export
##' @rdname ggstacy
#'
pggstacy <- function(q, alpha, gamma, kappa, log.p = FALSE, lower.tail = TRUE){
  x <- (q/gamma)^kappa
  p <- stats::pgamma(x, shape = alpha/gamma, scale = 1, lower.tail = lower.tail)
  if(isTRUE(log.p)){
      return(log(p))
    }else{
      return(p)
    }
}

##' @export
##' @rdname ggstacy
#'
qggstacy <- function(p, alpha = 1, gamma = 1, kappa = 1, log.p = FALSE, lower.tail = TRUE, ...){
  q <- stats::qgamma(p, shape = alpha/gamma, scale = 1, lower.tail = lower.tail)
  x <- gamma*(q^(1/kappa))
  if(isTRUE(log.p)){
    return(log(x))
  }else{
    return(x)
  }

}


##' @export
##' @rdname ggstacy
#'
rggstacy <- function(n, alpha = 1, gamma = 1, kappa = 1, ...){
  u <- stats::runif(n)
  x <- qggstacy(u, alpha, gamma, kappa)
  return(x)
}
