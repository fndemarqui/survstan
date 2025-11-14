# Probability function, distribution function, quantile function and random generation for the Piecewise Exponential (PE) distribution.

Probability function, distribution function, quantile function and
random generation for the Piecewise Exponential (PE) distribution.

## Usage

``` r
dpexp(x, rho, rates, log = FALSE)

ppexp(q, rho, rates, lower.tail = TRUE, log.p = FALSE)

qpexp(p, rho, rates, lower.tail = TRUE, log.p = FALSE)

rpexp(n, rho, rates)
```

## Arguments

- x:

  vector of time points.

- rho:

  vector of time grid knots.

- rates:

  vector of failure rates.

- log, log.p:

  logical; if TRUE, probabilities p are given as log(p).

- q:

  vector of quantiles.

- lower.tail:

  logical; if TRUE (default), probabilities are \\P\[X \le x\]\\;
  otherwise, \\P\[X \> x\]\\.

- p:

  vector of probabilities.

- n:

  number of random values to return.

## Value

dpexp gives the (log) probability function, ppexp gives the (log)
distribution function, qpexp gives the quantile function, and rpexp
generates random deviates.

## Examples

``` r
n <- 10
rho <- c(0, 1, 3, 7, Inf)
rates <- c(0.5, 4, 0.8, 0.1)
x <- sort(rpexp(n, rho=rho, rates=rates))
Fx <- ppexp(x, rho, rates)
y <- qpexp(Fx, rho, rates)
# checking:
x==y
#>  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
```
