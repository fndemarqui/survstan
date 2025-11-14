# Model.matrix method for survstan models

Reconstruct the model matrix for a survstan model.

## Usage

``` r
# S3 method for class 'survstan'
model.matrix(object, ...)
```

## Arguments

- object:

  an object of the class survstan.

- ...:

  further arguments passed to or from other methods.

## Value

The model matrix (or matrices) for the fit.

## Examples

``` r
# \donttest{
library(survstan)
fit <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
model.matrix(fit)
#>    ecog.ps rx
#> 1        1  1
#> 2        1  1
#> 3        2  1
#> 4        1  2
#> 5        1  1
#> 6        2  1
#> 7        2  2
#> 8        2  2
#> 9        1  1
#> 10       2  2
#> 11       2  1
#> 12       1  2
#> 13       2  2
#> 14       1  2
#> 15       1  1
#> 16       2  1
#> 17       2  1
#> 18       1  1
#> 19       1  2
#> 20       1  2
#> 21       2  2
#> 22       2  1
#> 23       1  1
#> 24       2  2
#> 25       1  2
#> 26       1  2
#> attr(,"assign")
#> [1] 1 2
# }
```
