# The Generalized Gamma Distribution (Prentice's alternative parametrization)

Probability function, distribution function, quantile function and
random generation for the distribution with parameters mu, sigma and
varphi.

## Usage

``` r
dggprentice(x, mu, sigma, varphi, log = FALSE)

pggprentice(q, mu = 0, sigma = 1, varphi, lower.tail = TRUE, log.p = FALSE)

qggprentice(p, mu = 0, sigma = 1, varphi, lower.tail = TRUE, log.p = FALSE)

rggprentice(n, mu = 0, sigma = 1, varphi, ...)
```

## Arguments

- x:

  vector of (non-negative integer) quantiles.

- mu:

  location parameter of the distribution.

- sigma:

  scale parameter of the distribution (sigma \> 0).

- varphi:

  shape parameter of the distribution.

- log, log.p:

  logical; if TRUE, probabilities p are given as log(p).

- q:

  vector of quantiles.

- lower.tail:

  logical; if TRUE (default), probabilities are \\P\[X \le x\]\\;
  otherwise, \\P\[X \> x\]\\.

- p:

  vector of probabilities.

- n:

  number of random values to return.

- ...:

  further arguments passed to other methods.

## Value

dggprentice gives the (log) probability function, pggprentice gives the
(log) distribution function, qggprentice gives the quantile function,
and rggprentice generates random deviates.

## Details

Probability density function: \$\$ f(x \| \mu, \sigma, \varphi) =
\begin{cases} \frac{\|\varphi\|(\varphi^{-2})^{\varphi^{-2}}}{\sigma
x\Gamma(\varphi^{-2})}\exp\\\varphi^{-2}\[\varphi w - \exp(\varphi
w)\]\\I\_{\[0, \infty)}(x), & \varphi \neq 0 \\
\frac{1}{\sqrt{2\pi}x\sigma}\exp\left\\-\frac{1}{2}\left(\frac{log(x)-\mu}{\sigma}\right)^2\right\\I\_{\[0,
\infty)}(x), & \varphi = 0 \end{cases} \$\$ where \\w = \frac{\log(x) -
\mu}{\sigma}\\, for \\-\infty \< \mu \< \infty\\, \\\sigma\>0\\ and
\\-\infty \< \varphi \< \infty\\.

Distribution function: \$\$ F(x\|\mu, \sigma, \varphi) = \begin{cases}
F\_{G}(y\|1/\varphi^2, 1), & \varphi \> 0 \\ 1-F\_{G}(y\|1/\varphi^2,
1), & \varphi \< 0 \\ F\_{LN}(x\|\mu, \sigma), & \varphi = 0 \end{cases}
\$\$ where \\y = \displaystyle\left(\frac{x}{\sigma}\right)^\varphi\\,
\\F\_{G}(\cdot\|\nu, 1)\\ is the distribution function of a gamma
distribution with shape parameter \\1/\varphi^2\\ and scale parameter
equals to 1, and \\F\_{LN}(x\|\mu, \sigma)\\ corresponds to the
distribution function of a lognormal distribution with location
parameter \\\mu\\ and scale parameter \\\sigma\\.
