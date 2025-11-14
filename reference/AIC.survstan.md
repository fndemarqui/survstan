# Akaike information criterion

Akaike information criterion

## Usage

``` r
# S3 method for class 'survstan'
AIC(object, ..., k = 2)
```

## Arguments

- object:

  an object of the class survstan.

- ...:

  further arguments passed to or from other methods.

- k:

  numeric, the penalty per parameter to be used; the default k = 2 is
  the classical AIC.

## Value

the Akaike information criterion value when a single model is passed to
the function; otherwise, a data.frame with the Akaike information
criterion values and the number of parameters is returned.

## Examples

``` r
# \donttest{
library(survstan)
fit1 <- aftreg(Surv(futime, fustat) ~ 1, data = ovarian, baseline = "weibull", init = 0)
fit2 <- aftreg(Surv(futime, fustat) ~ rx, data = ovarian, baseline = "weibull", init = 0)
fit3 <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
AIC(fit1, fit2, fit3)
#>    fit      aic npars
#> 1 fit1 199.9078     2
#> 2 fit2 200.7283     3
#> 3 fit3 202.1690     4
# }
```
