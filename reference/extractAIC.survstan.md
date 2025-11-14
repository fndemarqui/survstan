# Extract AIC from a Fitted Model

Computes the (generalized) Akaike An Information Criterion for a fitted
parametric model.

## Usage

``` r
# S3 method for class 'survstan'
extractAIC(fit, scale, k = 2, ...)
```

## Arguments

- fit:

  a fitted model of the class survstan

- scale:

  optional numeric specifying the scale parameter of the model.
  Currently only used in the "lm" method, where scale specifies the
  estimate of the error variance, and scale = 0 indicates that it is to
  be estimated by maximum likelihood.

- k:

  numeric specifying the ‘weight’ of the equivalent degrees of freedom
  part in the AIC formula.

- ...:

  further arguments passed to or from other methods.

## Value

the ANOVA table.

## Examples

``` r
# \donttest{
library(survstan)
fit1 <- aftreg(Surv(futime, fustat) ~ 1, data = ovarian, baseline = "weibull", init = 0)
fit2 <- aftreg(Surv(futime, fustat) ~ rx, data = ovarian, baseline = "weibull", init = 0)
fit3 <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
extractAIC(fit1)
#> [1]   2.0000 199.9078
extractAIC(fit2)
#> [1]   3.0000 200.7283
extractAIC(fit3)
#> [1]   4.000 202.169
# }
```
