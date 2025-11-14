# Variance-covariance matrix

This function extracts and returns the variance-covariance matrix
associated with the regression coefficients when the maximum likelihood
estimation approach is used in the model fitting.

## Usage

``` r
# S3 method for class 'survstan'
vcov(object, all = FALSE, ...)
```

## Arguments

- object:

  an object of the class survstan.

- all:

  logical; if FALSE (default), only covariance matrix associated with
  regression coefficients is returned; if TRUE, the full covariance
  matrix is returned.

- ...:

  further arguments passed to or from other methods.

## Value

the variance-covariance matrix associated with the parameters
estimators.

## Examples

``` r
# \donttest{
library(survstan)
fit <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
vcov(fit)
#>           ecog.ps        rx
#> ecog.ps 0.2777096 0.0110424
#> rx      0.0110424 0.2800500
# }
```
