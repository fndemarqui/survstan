# ggresiduals method for survstan models

ggresiduals method for survstan models

## Usage

``` r
# S3 method for class 'survstan'
ggresiduals(object, type = c("coxsnell", "martingale", "deviance"), ...)
```

## Arguments

- object:

  a fitted model object of the class survstan.

- type:

  type of residuals used in the plot: coxsnell (default), martingale and
  deviance.

- ...:

  further arguments passed to or from other methods.

## Value

the desired residual plot.

## Details

This function produces residuals plots of Cox-Snell residuals,
martingale residuals and deviance residuals.

## Examples

``` r
# \donttest{
library(survstan)
ovarian$rx <- as.factor(ovarian$rx)
fit <- aftreg(Surv(futime, fustat) ~ age + rx, data = ovarian, baseline = "weibull", init = 0)
ggresiduals(fit, type = "coxsnell")

ggresiduals(fit, type = "martingale")
#> `geom_smooth()` using method = 'loess' and formula = 'y ~ x'

ggresiduals(fit, type = "deviance")

# }
```
