# Fitting Proportional Odds Models

Function to fit proportional odds (PO) models.

## Usage

``` r
poreg(formula, data, baseline = "weibull", dist = NULL, init = 0, ...)
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

poreg returns an object of class "poreg" containing the fitted model.

## Examples

``` r
# \donttest{
library(survstan)
fit <- poreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull")
summary(fit)
#> Call:
#> poreg(formula = Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, 
#>     baseline = "weibull")
#> 
#> Proportional odds model fit with weibull baseline distribution: 
#> 
#> Regression coefficients:
#>         Estimate Std. Error z value Pr(>|z|)
#> ecog.ps  0.81190    0.71810  1.1306   0.2582
#> rx      -0.69296    0.71685 -0.9667   0.3337
#> 
#> Baseline parameters:
#>         Estimate Std. Error       2.5%     97.5%
#> alpha    1.18318    0.36436    0.64703    2.1636
#> gamma 1340.62000 1086.20339  273.93102 6561.0022
#> --- 
#> loglik = -97.03918   AIC = 202.0784 
# }

```
