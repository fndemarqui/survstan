# Fitting Proportional Hazards Models

Function to fit proportional hazards (PH) models.

## Usage

``` r
phreg(formula, data, baseline = "weibull", dist = NULL, init = 0, ...)
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

phreg returns an object of class "phreg" containing the fitted model.

## Examples

``` r
# \donttest{
library(survstan)
fit <- phreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull")
summary(fit)
#> Call:
#> phreg(formula = Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, 
#>     baseline = "weibull")
#> 
#> Proportional hazards model fit with weibull baseline distribution: 
#> 
#> Regression coefficients:
#>         Estimate Std. Error z value Pr(>|z|)
#> ecog.ps  0.43563    0.58727  0.7418   0.4582
#> rx      -0.59811    0.58823 -1.0168   0.3093
#> 
#> Baseline parameters:
#>         Estimate Std. Error       2.5%     97.5%
#> alpha    1.13138    0.28535    0.69012    1.8548
#> gamma  988.99454 1164.58729   98.37037 9943.1386
#> --- 
#> loglik = -97.08449   AIC = 202.169 
# }

```
