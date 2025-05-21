


#' Probability function, distribution function, quantile function and random generation for the Piecewise Exponential (PE) distribution.
#' @name pexp
#' @aliases pexp
#' @param x vector of time points.
#' @param q	 vector of quantiles.
#' @param p	 vector of probabilities.
#' @param n	 number of random values to return.
#' @param rho vector of time grid knots.
#' @param rates vector of failure rates.
#' @param log,log.p	 logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	logical; if TRUE (default), probabilities are \eqn{P[X \le x]}; otherwise, \eqn{P[X > x]}.
#' @return dpexp gives the (log) probability function, ppexp gives the (log) distribution function, qpexp gives the quantile function, and rpexp generates random deviates.
#' @description Probability function, distribution function, quantile function and random generation for the Piecewise Exponential (PE) distribution.
#'
#' @examples
#' n <- 10
#' rho <- c(0, 1, 3, 7, Inf)
#' rates <- c(0.5, 4, 0.8, 0.1)
#' x <- sort(rpexp(n, rho=rho, rates=rates))
#' Fx <- ppexp(x, rho, rates)
#' y <- qpexp(Fx, rho, rates)
#' # checking:
#' x==y

#' @rdname pexp
#' @export
dpexp <- function(x, rho, rates, log=FALSE){
  Ht <- Hpexp(x, rho, rates)
  id <- as.numeric(cut(x, rho, include.lowest = TRUE))
  log_ft <- log(rates[id]) - Ht
  if(log==FALSE)
    return(exp(log_ft))
  else
    return(log_ft)
}

#' @rdname pexp
#' @export
#'
ppexp <- function(q, rho, rates, lower.tail=TRUE, log.p=FALSE){
  Ht <- Hpexp(q, rho, rates)
  if(isTRUE(lower.tail)){
    p <- -expm1(-Ht)
  }else{
    p <- exp(-Ht)
  }
  if(isTRUE(log.p)){
    return(log(p))
  }else{
    return(p)
  }
}


#' @rdname pexp
#' @export
#'
#'
qpexp <- function(p, rho, rates, lower.tail=TRUE, log.p=FALSE){
  if(lower.tail==TRUE){
    u <- p
  }else{
    u <- 1-p
  }

  h <- diff(rho)
  n <- length(p)
  area <- c(0,cumsum(h*rates))
  Ft <- 1-exp(-area)
  id <- as.numeric(cut(u, unique(Ft), include.lowest = TRUE))

  if(n==1){
    if(id==1){
      q <- -log(1-u)/rates[1]
    }else{
      q <- rho[id] - (area[id] + log(1-u))/rates[id]
    }
  }else{
    q  <- rep(0, n)
    q[which(id==1)] <- -log(1-u[which(id==1)])/rates[1]
    aux <- which(id>1)
    q[aux] <- rho[id[aux]] - (area[id[aux]] + log(1-u[aux]))/rates[id[aux]]
  }
  return(q)
}


#' @rdname pexp
#' @export
#'
rpexp <- function(n, rho, rates){
  u <- stats::runif(n)
  x <- qpexp(u, rho, rates)
  return(x)
}


#---------------------------------------------

#' Hazard and cumulative hazard functions of the PE distribution
#' @name pehaz
#' @param x vector of time points.
#' @param rho vector of time grid knots.
#' @param rates vector of failure rates.
#' @return hpexp gives the hazard function and Hpexp gives the cumulative hazard function of the PE distribution.
#'

#' @rdname pehaz
#' @export
#'
hpexp <- function(x, rho, rates){
  id <- as.numeric(cut(x, rho, include.lowest = TRUE))
  return(rates[id])
}


#' @rdname pehaz
#' @export
#'
Hpexp <- function(x, rho, rates){
  h <- diff(rho)
  n <- length(x)
  id <- as.numeric(cut(x, rho, include.lowest = TRUE))
  area <- c(0, cumsum(h*rates))
  Ht <- (x-rho[id])*rates[id] + area[id]
  return(Ht)
}

#----------------------------------

#' Time grid
#' @export
#' @param time Vector of failure times
#' @param event Vector of failure indicators
#' @param m Optional. Number of intervals. If \code{NULL}, the number of intervals is
#' set to be equal to the number of distinct observed failure times.
#' @return Time grid.

time_grid <- function(time, event, m=NULL){
  n_int <- m
  o <- order(time)
  time <- time[o]
  event <- event[o]
  time.aux <- unique(time[event==1])
  if(is.null(n_int)){
    n_int <- length(time.aux)
  }

  m <- length(time.aux)
  if(n_int > m){
    rho <- c(0,unique(time[event==1]))
    rho[length(rho)] <- Inf
  }
  else{
    b <- min(m,n_int)
    k1 <- trunc(m/b)
    r <- m-b*k1
    k2 <- k1+1
    idf1 <- seq(k1,(b-r)*k1, k1)
    idf2 <- sort(seq(m,max(idf1),-k2))
    idf <- unique(c(idf1,idf2))
    rho <- c(0,time.aux[idf])
    rho[length(rho)] <- Inf
  }
  return(rho)
}


findInt <- function(time, rho){
  id <- as.numeric(cut(time, rho, include.lowest = TRUE))
  return(id)
}


#' Piecewise exponential baseline
#' @description This function is used to allow the user to specify the piecewise exponential baseline with arbitrary time grid/number of intervals.
#'
#' @aliases piecewise
#' @export
#' @param rho the specified time grid; default is NULL.
#' @param m the number of intervals; default is NULL.
#' @return a list with the baseline name, and the provided time grid and number of intervals.
piecewise <- function(rho = NULL, m = NULL){
  return(list(baseline = "piecewise", rho = rho, m = m))
}
