# Estimated standard errors

Estimated standard errors

## Usage

``` r
# S3 method for class 'survstan'
se(object, ...)
```

## Arguments

- object:

  an object of the class survstan.

- ...:

  further arguments passed to or from other methods.

## Value

a vector with the standard errors.

## Examples

``` r
# \donttest{
library(survstan)
fit <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
se(fit)
#>   ecog.ps        rx 
#> 0.5269816 0.5291975 
# }
```
