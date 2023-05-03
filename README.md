
<!-- README.md is generated from README.Rmd. Please edit that file -->

# survstan

<!-- badges: start -->
<!-- badges: end -->

The aim of the R package survstan is to provide a toolkit for fitting
survival models using Stan.

The R package survstan can be used to fit right-censored survival data
under independent censoring. The implemented models allow the fitting of
survival data in the presence/absence of covariates. All inferential
procedures are currently based on the maximum likelihood (ML) approach.

## Installation

You can install the released version of survstan from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("survstan")
```

You can install the development version of survstan from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("fndemarqui/survstan")
```

## Inference procedures

Let $(t_{i}, \delta_{i})$ be the observed survival time and its
corresponding failure indicator, $i=1, \cdots, n$, and
$\boldsymbol{\theta}$ be a $k \times 1$ vector of parameters. Then, the
likelihood function for right-censored survival data under independent
censoring can be expressed as:

The maximum likelihood estimate (MLE) of $\boldsymbol{\theta}$ is
obtained by directly maximization of $\log(L(\boldsymbol{\theta}))$
using the `rstan::optimizing()` function. The function
`rstan::optimizing()` further provides the hessian matrix of
$\log(L(\boldsymbol{\theta}))$, needed to obtain the observed Fisher
information matrix, which is given by:

Inferences on $\boldsymbol{\theta}$ are then based on the asymptotic
properties of the MLE, $\hat{\boldsymbol{\theta}}$, that state that:

## Baseline Distributions

Some of the most popular baseline survival distributions are implemented
in the R package survstan. Such distributions include:

- Exponential
- Weibull
- Lognormal
- Loglogistic

The parametrizations adopted in the package survstan are presented next.

### Exponential Distribution

If $T \sim \mbox{Exp}(\lambda)$, then

$$
        f(t|\lambda) = \lambda\exp\left\{-\lambda t\right\}I_{[0, \infty)}(t),
$$ where $\lambda>0$ is the rate parameter.

The survival and hazard functions in this case are given by:

$$
        S(t|\lambda) = \exp\left\{-\lambda t\right\}
$$ and $$
        h(t|\lambda) = \lambda.
$$

### Weibull Distribution

If $T \sim \mbox{Weibull}(\alpha, \gamma)$, then

$$
f(t|\alpha, \gamma) = \frac{\alpha}{\gamma^{\alpha}}t^{\alpha-1}\exp\left\{-\left(\frac{t}{\gamma}\right)^{\alpha}\right\}I_{[0, \infty)}(t),
$$ where $\alpha>0$ and $\gamma>0$ are the shape and scale parameters,
respectively.

The survival and hazard functions in this case are given by:

$$
        S(t|\alpha, \gamma) = \exp\left\{-\left(\frac{t}{\gamma}\right)^{\alpha}\right\}
$$ and $$
        h(t|\alpha, \gamma) = \frac{\alpha}{\gamma^{\alpha}}t^{\alpha-1}.
$$

### Lognormal Distribution

If $T \sim \mbox{LN}(\mu, \sigma)$, then

$$
f(t|\mu, \sigma) = \frac{1}{\sqrt{2\pi}t\sigma}\exp\left\{-\frac{1}{2}\left(\frac{log(t)-\mu}{\sigma}\right)^2\right\}I_{[0, \infty)}(t),
$$ where $-\infty < \mu < \infty$ and $\sigma>0$ are the mean and
standard deviation in the log scale of $T$.

The survival and hazard functions in this case are given by:

$$S(t|\mu, \sigma) = \Phi\left(\frac{-log(t)+\mu}{\sigma}\right)$$ and
$$h(t|\mu, \sigma) = \frac{f(t|\mu, \sigma)}{S(t|\mu, \sigma)},$$ where
$\Phi(\cdot)$ is the cumulative distribution function of the standard
normal distribution.

### Loglogistic Distribution

If $T \sim \mbox{LL}(\alpha, \gamma)$, then

$$
    f(t|\alpha, \gamma) = \frac{\frac{\alpha}{\gamma}\left(\frac{t}{\gamma}\right)^{\alpha-1}}{\left[1 + \left(\frac{t}{\gamma}\right)^{\alpha}\right]^2}I_{[0, \infty)}(t), ~ \alpha>0, \gamma>0,
$$

where $\alpha>0$ and $\gamma>0$ are the shape and scale parameters,
respectively.

The survival and hazard functions in this case are given by:

$$S(t|\alpha, \gamma) = \frac{1}{1+ \left(\frac{t}{\gamma}\right)^{\alpha}}$$
and $$
    h(t|\alpha, \gamma) = \frac{\frac{\alpha}{\gamma}\left(\frac{t}{\gamma}\right)^{\alpha-1}}{1 + \left(\frac{t}{\gamma}\right)^{\alpha}}.
$$

## Regression models

When covariates are available, it is possible to fit four different
regression models with the R package survstan:

- accelerated failure time (AFT) models;
- proportional hazards (PH) models;
- proportional odds (PO) models;
- accelerated hazard (AH) models.

Let $\mathbf{x}$ be a $1\times p$ vector of covariates,
$\boldsymbol{\beta}$ a $p \times 1$ of regression coefficients, and
$\boldsymbol{\theta}$ a vector of parameters associated with some
baseline survival distribution, and denote by
$\boldsymbol{\Theta} = (\boldsymbol{\theta}, \boldsymbol{\beta})^{T}$
the full vector of parameters. Here, to ensure identifiability, in all
regression structures the linear predictor
$\mathbf{x} \boldsymbol{\beta}$ does not include a intercept term.

The regression survival models implemented in the R package survstan are
briefly described in the sequel.

### Accelerate Failure Time Models

Accelerated failure time (AFT) models are defined as

$$
T = \exp\{\mathbf{x} \boldsymbol{\beta}\}\nu,
$$ where $\nu$ follows a baseline distribution with survival function
$S_{0}(\cdot|\boldsymbol{\theta})$ so that

$$
f(t|\boldsymbol{\Theta}, \mathbf{x}) = e^{-\mathbf{x} \boldsymbol{\beta}}f_{0}(te^{-\mathbf{x} \boldsymbol{\beta}}|\boldsymbol{\theta})
$$ and

$$
S(t|\boldsymbol{\Theta}, \mathbf{x}) = S_{0}(t e^{-\mathbf{x} \boldsymbol{\beta}}|\boldsymbol{\theta}).
$$

### Proportional Hazards Models

Proportional hazards (PH) models are defined as

$$
h(t|\Theta, \mathbf{x}) = h_{0}(t|\boldsymbol{\theta})\exp\{\mathbf{x} \boldsymbol{\beta}\},
$$ where $h_{0}(t|\boldsymbol{\theta})$ is a baseline hazard function so
that

$$
f(t|\boldsymbol{\Theta}, \mathbf{x}) = h_{0}(t|\boldsymbol{\theta})\exp\left\{\mathbf{x} \boldsymbol{\beta} - H_{0}(t|\boldsymbol{\theta})e^{\mathbf{x} \boldsymbol{\beta}}\right\},
$$ and

$$
S(t|\boldsymbol{\Theta}, \mathbf{x}) = \exp\left\{ - H_{0}(t|\boldsymbol{\theta})e^{\mathbf{x} \boldsymbol{\beta}}\right\}.
$$

### Proportional Odds Models

Proportional Odds (PO) models are defined as

$$
R(t|\Theta, \mathbf{x}) = R_{0}(t|\boldsymbol{\theta})\exp\{\mathbf{x} \boldsymbol{\beta}\},
$$ where
$\displaystyle R_{0}(t|\boldsymbol{\theta}) = \frac{1-S_{0}(t|\boldsymbol{\theta})}{S_{0}(t|\boldsymbol{\theta})} = \exp\{H_{0}(t|\boldsymbol{\theta})\}-1$
is a baseline odds function so that

$$
f(t|\boldsymbol{\Theta}, \mathbf{x}) = \frac{h_{0}(t|\boldsymbol{\theta})\exp\{\mathbf{x} \boldsymbol{\beta} + H_{0}(t|\boldsymbol{\theta})\}}{[1 + R_{0}(t|\boldsymbol{\theta})e^{\mathbf{x} \boldsymbol{\beta}}]^2}.
$$

and

$$
S(t|\boldsymbol{\Theta}, \mathbf{x}) = \frac{1}{1 + R_{0}(t|\boldsymbol{\theta})e^{\mathbf{x} \boldsymbol{\beta}}}.
$$

### Accelerated Hazard Models

Accelerated hazard (AH) models are defined as

$$h(t|\boldsymbol{\Theta},\mathbf{x}) = h_{0}\left(te^{\mathbf{x}\boldsymbol{\beta}}|\boldsymbol{\theta}\right)$$

so that

$$S(t|\boldsymbol{\Theta},\mathbf{x}) = \exp\left\{- H_{0}\left(t e^{\mathbf{x}\boldsymbol{\beta}}|\boldsymbol{\theta}\right)e^{-\mathbf{x}\boldsymbol{\beta}}
\right\}
$$ and
$$f(t|\boldsymbol{\theta}, \mathbf{x}) = h_{0}\left(te^{\mathbf{x}\boldsymbol{\beta}}|\boldsymbol{\theta}\right)\exp\left\{- H_{0}\left(t e^{\mathbf{x}\boldsymbol{\beta}}|\boldsymbol{\theta}\right)e^{-\mathbf{x}\boldsymbol{\beta}}
\right\}.
$$
