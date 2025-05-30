
# max allowed value of m when it is left unspecified by the user.
m_max <- 15


#' Bernstein polynomial
#' @description This function is used to allow the user to specify an arbitrary value for the polynomial's degree m. If m = NULL, then m = min(m_max, ceiling(n^0.4)) is used, where m_max = 15.
#'
#' @aliases bernstein
#' @export
#' @param m the Bernstein polynomial's degree; default is NULL.
#' @return a list with the baseline name and the polynomial's degree m.
bernstein <- function(m = NULL){
  return(list(baseline = "bernstein", m = m))
}


bp_b <- function(x, m, tau){
  n <- length(x)
  y <- x/tau
  b <- matrix(nrow=n, ncol=m)
  for(k in 1:m){
    b[,k] <- stats::dbeta(y, k, m - k + 1)/tau
  }
  return(b)
}

bp_B <- function(x, m, tau){
  n <- length(x)
  y <- x/tau
  B <- matrix(nrow=n, ncol=m)
  for(k in 1:m){
    B[,k] <- stats::pbeta(y, k, m - k + 1)
  }
  return(B)
}

# Computes the Bernstein polynomial's bases. (note: for computation stability, b is not divided by tau here)
BP <- function(x, m, tau) {
  #tau <- tau*(1)
  n <- length(x)
  y <- x/tau
  b <- matrix(nrow=n, ncol=m)
  B <- matrix(nrow=n, ncol=m)
  for(k in 1:m){
    b[,k] <- stats::dbeta(y, k, m - k + 1)
    B[,k] <- stats::pbeta(y, k, m - k + 1)
  }
  return(list(b=b, B=B))
}


dbernstein <- function(x, xi, tau, log = FALSE){
  #tau <- tau*(1+0.00001)
  m <- length(xi)
  bp <- BP(x, m, tau)
  lht <- as.numeric(log(bp$b%*%xi)) - log(tau)
  Ht <- as.numeric(bp$B%*%xi)
  lpdf <- lht - Ht
  if(log==FALSE){
    return(exp(lpdf))
  }else{
    return(lpdf)
  }
}


Hbernstein <- function(x, xi, tau){
  #tau <- tau*(1+0.00001)
  m <- length(xi)
  B <- bp_B(x, m, tau)
  Ht <- as.numeric(B%*%xi)
  return(Ht)
}

pbernstein <- function(x, xi, tau, lower.tail=TRUE, log.p=FALSE){
  #tau <- tau*(1+0.00001)
  m <- length(xi)
  B <- bp_B(x, m, tau)
  Ht <- as.numeric(B%*%xi)
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


qbernstein <- function(p, xi, tau, lower.tail = FALSE, log.p = FALSE, ...){

  #tau <- tau*(1+0.00001)

  if(isTRUE(lower.tail)){
    u <- 1-p
  }else{
    u <- p
  }

  cumhaz <- function(x, xi, tau){
    m <- length(xi)
    B <- bp_B(x, m, tau)
    Ht <- as.numeric(B%*%xi)
    return(Ht)
  }

  froot <- function(x, u, xi, tau){
    return( cumhaz(x, xi, tau) + log(u) )
  }
  q <- try(stats::uniroot(froot, interval = c(1.e-14, tau), tol = 1.e-4, u=u, xi=xi, tau=tau), TRUE)$root
  return(q)
}


rbernstein <- function(n, xi, tau){
  u <- stats::runif(n)
  x <- vector(length = n)
  for(i in 1:n){
    x[i] <- qbernstein(u[i], xi, tau)
  }
  return(x)
}

