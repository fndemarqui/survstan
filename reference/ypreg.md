# Fitting Yang and Prentice Models

Function to fit Yang and Prentice (YP) models.

## Usage

``` r
ypreg(formula, data, baseline = "weibull", dist = NULL, init = 0, ...)
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

ypreg returns an object of class "ypreg" containing the fitted model.

## Examples

``` r
# \donttest{
library(survstan)
fit <- ypreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull")
summary(fit)
#> Call:
#> ypreg(formula = Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, 
#>     baseline = "weibull")
#> 
#> Yang & Prentice model fit with weibull baseline distribution: 
#> 
#> Regression coefficients:
#>               Estimate Std. Error z value Pr(>|z|)
#> short-ecog.ps  1.36330    0.98921  1.3782   0.1681
#> short-rx      -1.21389    0.96329 -1.2601   0.2076
#> long-ecog.ps  -1.27129    1.40225 -0.9066   0.3646
#> long-rx        1.46122    2.35451  0.6206   0.5349
#> 
#> Baseline parameters:
#>         Estimate Std. Error       2.5%     97.5%
#> alpha    1.24544    0.47920    0.58588    2.6475
#> gamma 1265.53342 1026.11544  258.28865 6200.7171
#> --- 
#> loglik = -96.76574   AIC = 205.5315 
# }

```
