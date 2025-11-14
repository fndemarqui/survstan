# The Gompertz Distribution

Probability function, distribution function, quantile function and
random generation for the distribution with parameters alpha and gamma.

## Usage

``` r
dgompertz(x, alpha = 1, gamma = 1, log = FALSE, ...)

pgompertz(q, alpha = 1, gamma = 1, lower.tail = TRUE, log.p = FALSE, ...)

qgompertz(p, alpha = 1, gamma = 1, lower.tail = FALSE, log.p = FALSE, ...)

rgompertz(n, alpha = 1, gamma = 1, ...)
```

## Arguments

- x:

  vector of (non-negative integer) quantiles.

- alpha:

  shape parameter of the distribution (alpha \> 0).

- gamma:

  scale parameter of the distribution (gamma \> 0).

- log, log.p:

  logical; if TRUE, probabilities p are given as log(p).

- ...:

  further arguments passed to other methods.

- q:

  vector of quantiles.

- lower.tail:

  logical; if TRUE (default), probabilities are \\P\[X \le x\]\\;
  otherwise, \\P\[X \> x\]\\.

- p:

  vector of probabilities.

- n:

  number of random values to return.

## Value

dgompertz gives the (log) probability function, pgompertz gives the
(log) distribution function, qgompertz gives the quantile function, and
rgompertz generates random deviates.

## Details

Probability density function: \$\$ f(x\|\alpha, \gamma) = \alpha\gamma
\exp\\\gamma x - \alpha(e^{\gamma x} - 1)\\I\_{\[0, \infty)}(x), \$\$
for \\\alpha\>0\\ and \\\gamma\>0\\.

Distribution function: \$\$ F(x\|\alpha, \gamma) = 1 - \exp\\-
\alpha(e^{\gamma x} - 1)\\, \$\$ for \\x\>0\\, \\\alpha\>0\\ and
\\\gamma\>0\\.
