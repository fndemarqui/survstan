# Tidy a survstan object

Tidy a survstan object

## Usage

``` r
# S3 method for class 'survstan'
tidy(x, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE, ...)
```

## Arguments

- x:

  a fitted model object.

- conf.int:

  Logical indicating whether or not to include a confidence interval in
  the tidied output. Defaults to FALSE.

- conf.level:

  the confidence level required.

- exponentiate:

  Logical indicating whether or not to exponentiate the the coefficient
  estimates. Defaults to FALSE.

- ...:

  further arguments passed to or from other methods.

## Value

a tibble with a summary of the fit.

## Details

Convert a fitted model into a tibble.

## Examples

``` r
# \donttest{
library(survstan)
fit <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull")
tidy(fit)
#> # A tibble: 2 × 5
#>   term     estimate std.error statistic  p.value
#>   <chr>       <dbl>     <dbl>     <dbl>    <dbl>
#> 1 ecog.ps -0.385045  0.526982 -0.730662 0.464986
#> 2 rx       0.528758  0.529197  0.999169 0.317713
# }
```
