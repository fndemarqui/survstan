# Extract Log-Likelihood from a Fitted Model

Extracts the log-likelihood function for a fitted parametric model.

## Usage

``` r
# S3 method for class 'survstan'
logLik(object, ...)
```

## Arguments

- object:

  a fitted model of the class survstan

- ...:

  further arguments passed to or from other methods.

## Value

the log-likelihood value when a single model is passed to the function;
otherwise, a data.frame with the log-likelihood values and the number of
parameters is returned.

## Examples

``` r
# \donttest{
library(survstan)
fit1 <- aftreg(Surv(futime, fustat) ~ 1, data = ovarian, baseline = "weibull", init = 0)
fit2 <- aftreg(Surv(futime, fustat) ~ rx, data = ovarian, baseline = "weibull", init = 0)
fit3 <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
logLik(fit1, fit2, fit3)
#>    fit    loglik npars
#> 1 fit3 -97.08449     4
#> 2 fit2 -97.36415     3
#> 3 fit1 -97.95390     2
# }
```
