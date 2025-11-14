# anova method for survstan models

Compute analysis of variance (or deviance) tables for one or more fitted
model objects.

## Usage

``` r
# S3 method for class 'survstan'
anova(...)
```

## Arguments

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
anova(fit1, fit2, fit3)
#> 
#> weibull(aft) 1: Surv(futime, fustat) ~ 1 
#> weibull(aft) 2: Surv(futime, fustat) ~ rx 
#> weibull(aft) 3: Surv(futime, fustat) ~ ecog.ps + rx 
#> --- 
#>                    loglik        LR df Pr(>Chi)
#> weibull(aft) 1: -97.95390   1.73882  2   0.4192
#> weibull(aft) 2: -97.36415   0.55931  1   0.4545
#> weibull(aft) 3: -97.08449         -  -        -
# }
```
