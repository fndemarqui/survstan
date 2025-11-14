# Estimated regression coefficients

Estimated regression coefficients

## Usage

``` r
# S3 method for class 'survstan'
coef(object, ...)
```

## Arguments

- object:

  an object of the class survstan

- ...:

  further arguments passed to or from other methods

## Value

the estimated regression coefficients

## Examples

``` r
# \donttest{
library(survstan)
fit <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
coef(fit)
#>    ecog.ps         rx 
#> -0.3850453  0.5287577 
# }
```
