# The Generalized Gamma Distribution (Stacy's original parametrization)

Probability function, distribution function, quantile function and
random generation for the distribution with parameters alpha, gamma and
kappa.

## Usage

``` r
dggstacy(x, alpha, gamma, kappa, log = FALSE)

pggstacy(q, alpha, gamma, kappa, log.p = FALSE, lower.tail = TRUE)

qggstacy(
  p,
  alpha = 1,
  gamma = 1,
  kappa = 1,
  log.p = FALSE,
  lower.tail = TRUE,
  ...
)

rggstacy(n, alpha = 1, gamma = 1, kappa = 1, ...)
```

## Arguments

- x:

  vector of (non-negative integer) quantiles.

- alpha:

  shape parameter of the distribution (alpha \> 0).

- gamma:

  scale parameter of the distribution (gamma \> 0).

- kappa:

  shape parameter of the distribution (kappa \> 0).

- log, log.p:

  logical; if TRUE, probabilities p are given as log(p).

- q:

  vector of quantiles.

- lower.tail:

  logical; if TRUE (default), probabilities are \\P\[X \le x\]\\;
  otherwise, \\P\[X \> x\]\\.

- p:

  vector of probabilities.

- ...:

  further arguments passed to other methods.

- n:

  number of random values to return.

## Value

dggstacy gives the (log) probability function, pggstacy gives the (log)
distribution function, qggstacy gives the quantile function, and
rggstacy generates random deviates.

## Details

Probability density function: \$\$ f(x\|\alpha, \gamma, \kappa) =
\frac{\kappa}{\gamma^{\alpha}\Gamma(\alpha/\kappa)}x^{\alpha-1}\exp\left\\-\left(\frac{x}{\gamma}\right)^{\kappa}\right\\I\_{\[0,
\infty)}(x), \$\$ for \\\alpha\>0\\, \\\gamma\>0\\ and \\\kappa\>0\\.

Distribution function: \$\$ F(t\|\alpha, \gamma, \kappa) =
F\_{G}(x\|\nu, 1), \$\$ where \\x =
\displaystyle\left(\frac{t}{\gamma}\right)^\kappa\\, and
\\F\_{G}(\cdot\|\nu, 1)\\ corresponds to the distribution function of a
gamma distribution with shape parameter \\\nu = \alpha/\gamma\\ and
scale parameter equals to 1.
