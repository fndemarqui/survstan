# Support Functions for emmeans

Functions required for compatibility of survstan with emmeans. Users are
not required to call these functions themselves. Instead, they will be
called automatically by the `emmeans` function of the emmeans package.

## Usage

``` r
recover_data.survstan(object, ...)

recover_data.ypreg(object, term = c("short", "long"), ...)

recover_data.ehreg(object, term = c("AF", "RH"), ...)
```

## Arguments

- object:

  An object of the same class as is supported by a new method.

- ...:

  Additional parameters that may be supported by the method.

- term:

  character specifying whether AF or RH term regression coefficients are
  to be used.
