# residuals method for survstan models

residuals method for survstan models

## Usage

``` r
# S3 method for class 'survstan'
residuals(object, type = c("coxsnell", "martingale", "deviance"), ...)
```

## Arguments

- object:

  a fitted model object of the class survstan.

- type:

  type of residuals desired: coxsnell (default), martingale and
  deviance.

- ...:

  further arguments passed to or from other methods.

## Value

a vector containing the desired residuals.

## Details

This function extracts the residuals, martingale residuals and deviance
residuals of a survstan object.

## Examples

``` r
# \donttest{
library(survstan)
ovarian$rx <- as.factor(ovarian$rx)
fit <- aftreg(Surv(futime, fustat) ~ age + rx, data = ovarian, baseline = "weibull", init = 0)
residuals(fit, type = "coxsnell")
#>          1          2          3          4          5          6          7 
#> 0.10688361 0.48993205 0.26908236 0.08873277 0.16806224 0.43241247 0.17686609 
#>          8          9         10         11         12         13         14 
#> 0.28063085 1.47404617 0.19518977 0.86115946 0.15636767 0.65177670 0.45113772 
#>         15         16         17         18         19         20         21 
#> 0.10611859 0.20681060 0.16076382 0.40804195 0.57536345 0.16099606 1.51372713 
#>         22         23         24         25         26 
#> 2.28089328 0.03657446 0.26538744 0.33530877 0.14774074 
residuals(fit, type = "martingale")
#>           1           2           3           4           5           6 
#>  0.89311639  0.51006795  0.73091764 -0.08873277  0.83193776 -0.43241247 
#>           7           8           9          10          11          12 
#>  0.82313391  0.71936915 -1.47404617  0.80481023  0.13884054 -0.15636767 
#>          13          14          15          16          17          18 
#> -0.65177670 -0.45113772 -0.10611859 -0.20681060 -0.16076382 -0.40804195 
#>          19          20          21          22          23          24 
#> -0.57536345 -0.16099606 -1.51372713 -1.28089328  0.96342554  0.73461256 
#>          25          26 
#>  0.66469123 -0.14774074 
residuals(fit, type = "deviance")
#>          1          2          3          4          5          6          7 
#>  1.6388401  0.6378411  1.0787216 -0.4212666  1.3794804 -0.9299596  1.3485018 
#>          8          9         10         11         12         13         14 
#>  1.0500914 -1.7170010  1.2876123  0.1458427 -0.5592275 -1.1417326 -0.9498818 
#>         15         16         17         18         19         20         21 
#> -0.4606921 -0.6431339 -0.5670341 -0.9033736 -1.0727194 -0.5674435 -1.7399581 
#>         22         23         24         25         26 
#> -0.9553283  2.1656313  1.0880735  0.9252159 -0.5435821 
# }
```
