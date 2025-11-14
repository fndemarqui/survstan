# Bernstein polynomial

This function is used to allow the user to specify an arbitrary value
for the polynomial's degree m. If m = NULL, then m = min(m_max,
ceiling(n^0.4)) is used, where m_max = 15.

## Usage

``` r
bernstein(m = NULL)
```

## Arguments

- m:

  the Bernstein polynomial's degree; default is NULL.

## Value

a list with the baseline name and the polynomial's degree m.
