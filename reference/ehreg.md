# Fitting Extended Hazard Models

Function to fit Extended Hazard (EH) models.

## Usage

``` r
ehreg(formula, data, baseline = "weibull", dist = NULL, init = 0, ...)
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

ehreg returns an object of class "ehreg" containing the fitted model.

## Examples

``` r
# \donttest{
library(survstan)
fit <- ehreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull")
summary(fit)
#> Call:
#> ehreg(formula = Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, 
#>     baseline = "weibull")
#> 
#> Extended hazard model fit with weibull baseline distribution: 
#> 
#> Regression coefficients:
#>             Estimate Std. Error z value Pr(>|z|)
#> AF-ecog.ps -0.061010   0.075973 -0.8030   0.4219
#> AF-rx      -0.017487   0.076093 -0.2298   0.8182
#> RH-ecog.ps  0.426947   0.577325  0.7395   0.4596
#> RH-rx      -0.600648   0.578552 -1.0382   0.2992
#> 
#> Baseline parameters:
#>         Estimate Std. Error       2.5%     97.5%
#> alpha    1.13146    0.28536    0.69018    1.8549
#> gamma  987.77515 1162.96582   98.28527 9927.2229
#> --- 
#> loglik = -97.08449   AIC = 206.169 
# }

```
