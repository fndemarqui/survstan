# Computes the crossing survival times

Computes the crossing survival times

## Usage

``` r
# S3 method for class 'survstan'
cross_time(
  object,
  newdata1,
  newdata2,
  conf.level = 0.95,
  nboot = 1000,
  cores = 1,
  ...
)
```

## Arguments

- object:

  an object of class survstan

- newdata1:

  a data frame containing the first set of explanatory variables

- newdata2:

  a data frame containing the second set of explanatory variables

- conf.level:

  level of the confidence/credible intervals

- nboot:

  number of bootstrap samples (default nboot=1000).

- cores:

  number of cores to be used in the bootstrap sampling; default is 1
  core;

- ...:

  further arguments passed to or from other methods.

## Value

the crossing survival time

## Examples

``` r
# \donttest{
library(survstan)
data(ipass)
fit <- ypreg(Surv(time, status)~arm, data=ipass, baseline = "weibull")
summary(fit)
#> Call:
#> ypreg(formula = Surv(time, status) ~ arm, data = ipass, baseline = "weibull")
#> 
#> Yang & Prentice model fit with weibull baseline distribution: 
#> 
#> Regression coefficients:
#>            Estimate Std. Error  z value  Pr(>|z|)    
#> short-arm  1.361063   0.182368   7.4633  8.44e-14 ***
#> long-arm  -1.365391   0.082943 -16.4619 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Baseline parameters:
#>       Estimate Std. Error     2.5%  97.5%
#> alpha 1.846062   0.060577 1.731071 1.9687
#> gamma 6.961753   0.166121 6.643658 7.2951
#> --- 
#> loglik = -2772.375   AIC = 5552.751 
newdata1 <- data.frame(arm=0)
newdata2 <- data.frame(arm=1)
tcross <- cross_time(fit, newdata1, newdata2, nboot = 10)
#> Please, be patient!!!
#> Bootstrap samples draw using 1 core
tcross
#>   Estimate Std. Error     2.5%    97.5%
#> 1 5.721489  0.3370274 5.001726 6.023493
# }
```
