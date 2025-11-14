# Rank a collection of survstan models

Rank a collection of survstan models

## Usage

``` r
rank_models(formula, data, survreg, baseline, dist = NULL, ...)
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

- survreg:

  survival regression models to be fitted (AFT, AH, PH, PO, YP and EH).

- baseline:

  baseline distributions to be fitted; options currently available are:
  exponential, weibull, lognormal, loglogistic and Birnbaum-Saunders
  (fatigue) distributions.

- dist:

  alternative way to specify the baseline distributions (for compability
  with the [`survreg`](https://rdrr.io/pkg/survival/man/survreg.html)
  function); default is NULL.

- ...:

  further arguments passed to other methods.

## Value

a tibble containing the fitted models ranked according to their AICs.

## Examples

``` r
# \donttest{
library(survstan)
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union

veteran <- veteran %>%
  mutate(across(c(trt, prior, celltype), as.factor))
fits <- rank_models(
  formula = Surv(time, status) ~ celltype+karno,
  data = veteran,
  survreg = c("aftreg", "ahreg", "phreg", "poreg", "ypreg", "ehreg"),
  baseline = c("exponential", "weibull", "lognormal", "loglogistic")
)
#> Warning: There were 2 warnings in `dplyr::mutate()`.
#> The first warning was:
#> ℹ In argument: `fit = purrr::pmap(...)`.
#> Caused by warning in `ahreg()`:
#> ! The AH model with baseline exponential distribution is non-identifiable! Please, choose another baseline distribution.
#> ℹ Run `dplyr::last_dplyr_warnings()` to see the 1 remaining warning.
# }

```
