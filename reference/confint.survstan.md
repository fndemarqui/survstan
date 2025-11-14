# Confidence intervals for the regression coefficients

Confidence intervals for the regression coefficients

## Usage

``` r
# S3 method for class 'survstan'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  an object of the class survstan.

- parm:

  a specification of which parameters are to be given confidence
  intervals, either a vector of numbers or a vector of names. If
  missing, all parameters are considered.

- level:

  the confidence level required.

- ...:

  further arguments passed to or from other methods.

## Value

100(1-alpha) confidence intervals for the regression coefficients.

## Examples

``` r
# \donttest{
library(survstan)
fit <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
confint(fit)
#>               2.5%        97.5%
#> ecog.ps -1.4179102    0.6478196
#> rx      -0.5084503    1.5659657
#> alpha    0.6901645    1.8548609
#> gamma   98.3660009 9940.3556002
# }
```
