
#' The Generalized Gamma Distribution (Prentice's alternative parametrization)
#' @name ggprentice
#' @aliases ggprentice dggprentice pggprentice qggprentice rggprentice
#'
#' @keywords distribution
#' @concept Univariate
#' @concept Continuous
#' @description Probability function, distribution function, quantile function and random generation for the  distribution with parameters mu, sigma and varphi.
#'
#' @param x	vector of (non-negative integer) quantiles.
#' @param q	vector of quantiles.
#' @param p	vector of probabilities.
#' @param n	number of random values to return.
#' @param mu location parameter of the distribution.
#' @param sigma scale parameter of the distribution (sigma > 0).
#' @param varphi shape parameter of the distribution.
#' @param log,log.p	 logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	logical; if TRUE (default), probabilities are \eqn{P[X \le x]}; otherwise, \eqn{P[X > x]}.
#' @param ... further arguments passed to other methods.
#'
#' @details
#'
#' Probability density function:
#' \deqn{
#' f(x | \mu, \sigma, \varphi) =
#' \begin{cases}
#' \frac{|\varphi|(\varphi^{-2})^{\varphi^{-2}}}{\sigma x\Gamma(\varphi^{-2})}\exp\{\varphi^{-2}[\varphi w - \exp(\varphi w)]\}I_{[0, \infty)}(x), & \varphi \neq 0 \\
#' \frac{1}{\sqrt{2\pi}x\sigma}\exp\left\{-\frac{1}{2}\left(\frac{log(x)-\mu}{\sigma}\right)^2\right\}I_{[0, \infty)}(x), & \varphi = 0
#' \end{cases}
#'  }
#' where \eqn{w = \frac{\log(x) - \mu}{\sigma}}, for \eqn{-\infty < \mu < \infty}, \eqn{\sigma>0} and \eqn{-\infty < \varphi < \infty}.
#'
#'
#' Distribution function:
#' \deqn{
#' F(x|\mu, \sigma, \varphi) =
#'  \begin{cases}
#' F_{G}(y|1/\varphi^2, 1), & \varphi > 0 \\
#' 1-F_{G}(y|1/\varphi^2, 1), & \varphi < 0 \\
#' F_{LN}(x|\mu, \sigma), & \varphi = 0
#' \end{cases}
#' }
#' where \eqn{y = \displaystyle\left(\frac{x}{\sigma}\right)^\varphi},
#' \eqn{F_{G}(\cdot|\nu, 1)} is the distribution function of
#' a gamma distribution with shape parameter \eqn{1/\varphi^2} and scale
#' parameter equals to 1, and \eqn{F_{LN}(x|\mu, \sigma)} corresponds to the
#' distribution function of a lognormal distribution with location parameter
#' \eqn{\mu} and scale parameter \eqn{\sigma}.
#' @return dggprentice gives the (log) probability function, pggprentice gives the (log) distribution function, qggprentice gives the quantile function, and rggprentice generates random deviates.
#'


##' @export
##' @rdname ggprentice
#'
dggprentice <- function(x, mu, sigma, varphi, log = FALSE){
  if(varphi == 0){
    lpdf = stats::dlnorm(x, meanlog = mu, sdlog = sigma)
  }else{
    w = (log(x) - mu) / sigma
    abs_q = abs(varphi)
    q2 = varphi^(-2)
    qw = varphi*w
    lpdf = - log(sigma) - log(x) + (1 - 2 * q2)*log(abs_q) + q2 * (qw - exp(qw)) - lgamma(q2)
  }

  if(isTRUE(log)){
    return(exp(lpdf))
  }else{
    return(lpdf)
  }
}

#' @export
#' @rdname ggprentice
#'
pggprentice <- function (q, mu = 0, sigma = 1, varphi, lower.tail = TRUE, log.p = FALSE) {
  if(varphi == 0){
    p <- stats::plnorm(q, meanlog = mu, sdlog = sigma)
  }else{
    w = (log(q) - mu)/sigma
    q2 = varphi^(-2)
    aux = q2*exp(varphi*w)
    if(varphi > 0){
      p = stats::pgamma(aux, q2, 1, lower.tail = FALSE)
    }else if(varphi < 0){
      p = stats::pgamma(aux, q2, 1)
    }
  }

  if(!isTRUE(lower.tail)){
    p <- 1 - p
  }
  if (isTRUE(log.p)) {
    return(log(p))
  }
  else {
    return(p)
  }
}

#' @export
#' @rdname ggprentice
#'
qggprentice <- function (p, mu = 0, sigma = 1, varphi, lower.tail = TRUE, log.p = FALSE){
  if(varphi != 0){
    q2 <- varphi^(-2)
    if(varphi > 0){
      x = stats::qgamma(p, q2, 1)
    }else{
      x = stats::qgamma(p, q2, 1, lower.tail = FALSE)
    }
    q <- mu + sigma*(log(x) + 2*log(abs(varphi)))/varphi
  }else{
    q <- log(stats::qlnorm(p, meanlog = mu, sdlog = sigma, lower.tail = lower.tail))
  }
  if (isTRUE(log.p)) {
    return(q)
  }
  else {
    return(exp(q))
  }
}


#' @export
#' @rdname ggprentice
#'
rggprentice <- function(n, mu = 0, sigma = 1, varphi, ...){
  u <- stats::runif(n)
  x <- qggprentice(u, mu, sigma, varphi)
  return(x)
}
