# Parameters estimates of a survstan model

Parameters estimates of a survstan model

## Usage

``` r
estimates(object, ...)
```

## Arguments

- object:

  an object of the class survstan.

- ...:

  further arguments passed to or from other methods.

## Value

the parameters estimates of a given survstan model.

## Examples

``` r
# \donttest{
library(survstan)
fit <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
estimates(fit)
#>     ecog.ps          rx       alpha       gamma 
#>  -0.3850453   0.5287577   1.1314412 988.8341760 
# }
```
