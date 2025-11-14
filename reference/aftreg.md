# Fitting Accelerated Failure Time Models

Function to fit accelerated failure time (AFT) models.

## Usage

``` r
aftreg(formula, data, baseline = "weibull", dist = NULL, init = 0, ...)
```

## Arguments

- formula:

  an object of class "formula" (or one that can be coerced to that
  class): a symbolic description of the model to be fitted.

- data:

  data an optional data frame, list or environment (or object coercible
  by as.data.frame to a data frame) containing the variables in the
  model. If not found in data, the variables are taken from
  environment(formula), typically the environment from which function is
  called.

- baseline:

  the chosen baseline distribution; options currently available are:
  exponential, weibull, lognormal, loglogistic and Birnbaum-Saunders
  (fatigue) distributions.

- dist:

  alternative way to specify the baseline distribution (for
  compatibility with the
  [`survreg`](https://rdrr.io/pkg/survival/man/survreg.html) function);
  default is NULL.

- init:

  initial values specification (default value is 0); see the detailed
  documentation for `init` in
  [`optimizing`](https://mc-stan.org/rstan/reference/stanmodel-method-optimizing.html).

- ...:

  further arguments passed to other methods.

## Value

aftreg returns an object of class "aftreg" containing the fitted model.

## Examples

``` r
# \donttest{
library(survstan)
fit <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull")
summary(fit)
#> Call:
#> aftreg(formula = Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, 
#>     baseline = "weibull")
#> 
#> Accelerated failure time model fit with weibull baseline distribution: 
#> 
#> Regression coefficients:
#>         Estimate Std. Error z value Pr(>|z|)
#> ecog.ps -0.38505    0.52698 -0.7307   0.4650
#> rx       0.52876    0.52920  0.9992   0.3177
#> 
#> Baseline parameters:
#>         Estimate Std. Error       2.5%     97.5%
#> alpha    1.13144    0.28536    0.69016    1.8549
#> gamma  988.83418 1164.33904   98.36600 9940.3556
#> --- 
#> loglik = -97.08449   AIC = 202.169 
# }
```
