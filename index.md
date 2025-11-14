# survstan

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

Let $\left( t_{i},\delta_{i} \right)$ be the observed survival time and
its corresponding failure indicator, $i = 1,\cdots,n$, and $\mathbf{θ}$
be a $k \times 1$ vector of parameters. Then, the likelihood function
for right-censored survival data under independent censoring can be
expressed as:

$$L({\mathbf{θ}}) = \prod\limits_{i = 1}^{n}f\left( t_{i}|{\mathbf{θ}} \right)^{\delta_{i}}S\left( t_{i}|{\mathbf{θ}} \right)^{1 - \delta_{i}}.$$

The maximum likelihood estimate (MLE) of $\mathbf{θ}$ is obtained by
directly maximization of $\log\left( L({\mathbf{θ}}) \right)$ using the
[`rstan::optimizing()`](https://mc-stan.org/rstan/reference/stanmodel-method-optimizing.html)
function. The function
[`rstan::optimizing()`](https://mc-stan.org/rstan/reference/stanmodel-method-optimizing.html)
further provides the hessian matrix of
$\log\left( L({\mathbf{θ}}) \right)$, needed to obtain the observed
Fisher information matrix, which is given by:

$$\mathcal{I}\left( \widehat{\mathbf{θ}} \right) = - \frac{\partial^{2}}{\partial{\mathbf{θ}}{\mathbf{θ}}\prime}\log L({\mathbf{θ}}) \mid_{{\mathbf{θ}} = \widehat{\mathbf{θ}}},$$

Inferences on $\mathbf{θ}$ are then based on the asymptotic properties
of the MLE, $\widehat{\mathbf{θ}}$, that state that:

$$\widehat{\mathbf{θ}} \asymp N_{k}\left( {\mathbf{θ}},\mathcal{I}^{- 1}\left( \widehat{\mathbf{θ}} \right) \right).$$

## Baseline Distributions

Some of the most popular baseline survival distributions are implemented
in the R package survstan. Such distributions include:

- Exponential
- Weibull
- Lognormal
- Loglogistic
- Gamma,
- Generalized Gamma (original Stacy’s parametrization)
- Generalized Gamma (alternative Prentice’s parametrization)
- Gompertz
- Rayleigh
- Birnbaum-Saunders (fatigue)

The parametrizations adopted in the package survstan are presented next.

### Exponential Distribution

If $T \sim \text{Exp}(\lambda)$, then

$$f\left( t|\lambda \right) = \lambda\exp\left\{ - \lambda t \right\} I_{\lbrack 0,\infty)}(t),$$
where $\lambda > 0$ is the rate parameter.

The survival and hazard functions in this case are given by:

$$S\left( t|\lambda \right) = \exp\left\{ - \lambda t \right\}$$ and
$$h\left( t|\lambda \right) = \lambda.$$

### Weibull Distribution

If $T \sim \text{Weibull}(\alpha,\gamma)$, then

$$f\left( t|\alpha,\gamma \right) = \frac{\alpha}{\gamma^{\alpha}}t^{\alpha - 1}\exp\left\{ - \left( \frac{t}{\gamma} \right)^{\alpha} \right\} I_{\lbrack 0,\infty)}(t),$$
where $\alpha > 0$ and $\gamma > 0$ are the shape and scale parameters,
respectively.

The survival and hazard functions in this case are given by:

$$S\left( t|\alpha,\gamma \right) = \exp\left\{ - \left( \frac{t}{\gamma} \right)^{\alpha} \right\}$$
and
$$h\left( t|\alpha,\gamma \right) = \frac{\alpha}{\gamma^{\alpha}}t^{\alpha - 1}.$$

### Lognormal Distribution

If $T \sim \text{LN}(\mu,\sigma)$, then

$$f\left( t|\mu,\sigma \right) = \frac{1}{\sqrt{2\pi}t\sigma}\exp\left\{ - \frac{1}{2}\left( \frac{log(t) - \mu}{\sigma} \right)^{2} \right\} I_{\lbrack 0,\infty)}(t),$$
where $- \infty < \mu < \infty$ and $\sigma > 0$ are the mean and
standard deviation in the log scale of $T$.

The survival and hazard functions in this case are given by:

$$S\left( t|\mu,\sigma \right) = \Phi\left( \frac{- log(t) + \mu}{\sigma} \right)$$
and
$$h\left( t|\mu,\sigma \right) = \frac{f\left( t|\mu,\sigma \right)}{S\left( t|\mu,\sigma \right)},$$
where $\Phi( \cdot )$ is the cumulative distribution function of the
standard normal distribution.

### Loglogistic Distribution

If $T \sim \text{LL}(\alpha,\gamma)$, then

$$f\left( t|\alpha,\gamma \right) = \frac{\frac{\alpha}{\gamma}\left( \frac{t}{\gamma} \right)^{\alpha - 1}}{\left\lbrack 1 + \left( \frac{t}{\gamma} \right)^{\alpha} \right\rbrack^{2}}I_{\lbrack 0,\infty)}(t),\ \alpha > 0,\gamma > 0,$$

where $\alpha > 0$ and $\gamma > 0$ are the shape and scale parameters,
respectively.

The survival and hazard functions in this case are given by:

$$S\left( t|\alpha,\gamma \right) = \frac{1}{1 + \left( \frac{t}{\gamma} \right)^{\alpha}}$$
and
$$h\left( t|\alpha,\gamma \right) = \frac{\frac{\alpha}{\gamma}\left( \frac{t}{\gamma} \right)^{\alpha - 1}}{1 + \left( \frac{t}{\gamma} \right)^{\alpha}}.$$

### Gamma Distribution

If $T \sim \text{Gamma}(\alpha,\lambda)$, then

$$f\left( t|\alpha,\lambda \right) = \frac{\lambda^{\alpha}}{\Gamma(\alpha)}t^{\alpha - 1}\exp\left\{ - \lambda t \right\} I_{\lbrack 0,\infty)}(t),$$

where $\Gamma(\alpha) = \int_{0}^{\infty}u^{\alpha - 1}\exp\{ - u\} du$
is the gamma function.

The survival function is given by

$$S\left( t|\alpha,\lambda \right) = 1 - \frac{\gamma^{*}(\alpha,\lambda t)}{\Gamma(\alpha)},$$
where $\gamma^{*}(\alpha,\lambda t)$ is the lower incomplete gamma
function, which is available only numerically. Finally, the hazard
function is expressed as:

$$h\left( t|\alpha,\lambda \right) = \frac{f\left( t|\alpha,\lambda \right)}{S\left( t|\alpha,\lambda \right)}.$$

### Generalized Gamma Distribution (original Stacy’s parametrization)

If $T \sim \text{ggstacy}(\alpha,\gamma,\kappa)$, then

$$f\left( t|\alpha,\gamma,\kappa \right) = \frac{\kappa}{\gamma^{\alpha}\Gamma(\alpha/\kappa)}t^{\alpha - 1}\exp\left\{ - \left( \frac{t}{\gamma} \right)^{\kappa} \right\} I_{\lbrack 0,\infty)}(t),$$
for $\alpha > 0$, $\gamma > 0$ and $\kappa > 0$.

It can be show that the survival function can be expressed as:

$$S\left( t|\alpha,\gamma,\kappa \right) = S_{G}\left( x|\nu,1 \right),$$
where $x = \left( \frac{t}{\gamma} \right)^{\kappa}$, and
$F_{G}\left( \cdot |\nu,1 \right)$ corresponds to the distribution
function of a gamma distribution with shape parameter
$\nu = \alpha/\gamma$ and scale parameter equals to 1.

Finally, the hazard function is expressed as:

$$h\left( t|\alpha,\gamma,\kappa \right) = \frac{f\left( t|\alpha,\gamma,\kappa \right)}{S\left( t|\alpha,\gamma,\kappa \right)}.$$

### Generalized Gamma Distribution (alternative Prentice’s parametrization)

If $T \sim \text{ggprentice}(\mu,\sigma,\varphi)$, then

$$f\left( t|\mu,\sigma,\varphi \right) = \begin{cases}
{\frac{|\varphi|\left( \varphi^{- 2} \right)^{\varphi^{- 2}}}{\sigma t\Gamma\left( \varphi^{- 2} \right)}\exp\{\varphi^{- 2}\left\lbrack \varphi w - \exp(\varphi w) \right\rbrack\} I_{\lbrack 0,\infty)}(t),} & {\varphi \neq 0} \\
{\frac{1}{\sqrt{2\pi}t\sigma}\exp\left\{ - \frac{1}{2}\left( \frac{log(t) - \mu}{\sigma} \right)^{2} \right\} I_{\lbrack 0,\infty)}(t),} & {\varphi = 0}
\end{cases}$$ where $w = \frac{\log(t) - \mu}{\sigma}$, for
$- \infty < \mu < \infty$, $\sigma > 0$ and
$- \infty < \varphi < \infty$\$.

It can be show that the survival function can be expressed as:

$$S\left( t|\mu,\sigma,\varphi \right) = \begin{cases}
{S_{G}\left( x|1/\varphi^{2},1 \right),} & {\varphi > 0} \\
{1 - S_{G}\left( x|1/\varphi^{2},1 \right),} & {\varphi < 0} \\
{S_{LN}\left( x|\mu,\sigma \right),} & {\varphi = 0}
\end{cases}$$ where $x = \frac{1}{\varphi^{2}}\exp\{\varphi w\}$,
$S_{G}\left( \cdot |1/\varphi^{2},1 \right)$ is the distribution
function of a gamma distribution with shape parameter $1/\varphi^{2}$
and scale parameter equals to 1, and $S_{LN}\left( x|\mu,\sigma \right)$
corresponds to the survival function of a lognormal distribution with
location parameter $\mu$ and scale parameter $\sigma$.

Finally, the hazard function is expressed as:

$$h\left( t|\alpha,\gamma,\kappa \right) = \frac{f\left( t|\alpha,\gamma,\kappa \right)}{S\left( t|\alpha,\gamma,\kappa \right)}.$$

### Gompertz Distribution

If $T \sim \text{Gamma}(\alpha,\gamma)$, then

$$f\left( t|\alpha,\lambda \right) = \alpha\exp\left\{ \gamma t - \frac{\alpha}{\gamma}\left( e^{\gamma t} - 1 \right) \right\} I_{\lbrack 0,\infty)}(t).$$

The survival and hazard functions are given, respectively, by

$$S\left( t|\alpha,\lambda \right) = \exp\left\{ - \frac{\alpha}{\gamma}\left( e^{\gamma t} - 1 \right) \right\}.$$
and

\$\$h(t\|\alpha, \lambda) = \alpha\exp\\\gamma t}.\$\$

### Rayleigh Distribution

Let $T \sim \text{rayleigh}(\sigma)$, where $\sigma > 0$ is a scale
parameter. Then, the density, survival and hazard functions are
respectively given by:

$$f\left( t|\sigma \right) = \frac{x}{\sigma^{2}}\exp\left\{ - \frac{x^{2}}{2\sigma^{2}} \right\},$$$$S\left( t|\sigma \right) = \exp\left\{ - \frac{x^{2}}{2\sigma^{2}} \right\}$$
and

$$h\left( t|\sigma \right) = \frac{x}{\sigma^{2}}.$$

### Birnbaum-Saunders (fatigue) Distribution

If $T \sim \text{fatigue}(\alpha,\gamma)$, then

$$f\left( t|\alpha,\gamma \right) = \frac{\sqrt{\frac{t}{\gamma}} + \sqrt{\frac{\gamma}{t}}}{2\alpha t}\phi\left( \sqrt{\frac{t}{\gamma}} + \sqrt{\frac{\gamma}{t}} \right)(t),\ \alpha > 0,\gamma > 0,$$

where $\phi( \cdot )$ is the probability density function of a standard
normal distribution, $\alpha > 0$ and $\gamma > 0$ are the shape and
scale parameters, respectively.

The survival function in this case is given by:

$$S\left( t|\alpha,\gamma \right) = \Phi\left( \sqrt{\frac{t}{\gamma}} - \sqrt{\frac{\gamma}{t}} \right)(t)$$,

where $\Phi( \cdot )$ is the cumulative distribution function of a
standard normal distribution. The hazard function is given by
$$h\left( t|\mu,\sigma \right) = \frac{f\left( t|\alpha,\gamma \right)}{S\left( t|\alpha,\gamma \right)}.$$

## Regression models

When covariates are available, it is possible to fit six different
regression models with the R package survstan:

- accelerated failure time (AFT) models;
- proportional hazards (PH) models;
- proportional odds (PO) models;
- accelerated hazard (AH) models.
- Yang and Prentice (YP) models.
- extended hazard (EH) models.

The regression survival models implemented in the R package survstan are
briefly described in the sequel. Denote by $\mathbf{x}$ a $1 \times p$
vector of covariates, and let $\mathbf{β}$ and $\mathbf{ϕ}$ be
$p \times 1$ vectors of regression coefficients, and $\mathbf{θ}$ a
vector of parameters associated with some baseline survival
distribution. To prevent identifiability issues, it is assumed that the
linear predictors $\mathbf{x}{\mathbf{β}}$ and $\mathbf{x}{\mathbf{ϕ}}$
do not include an intercept term.

### Accelerate Failure Time Models

Accelerated failure time (AFT) models are defined as

$$T = \exp\{\mathbf{x}{\mathbf{β}}\}\nu,$$ where $\nu$ follows a
baseline distribution with survival function
$S_{0}\left( \cdot |{\mathbf{θ}} \right)$ so that

$$f\left( t|{\mathbf{θ}},{\mathbf{β}},\mathbf{x} \right) = e^{- \mathbf{x}{\mathbf{β}}}f_{0}\left( te^{- \mathbf{x}{\mathbf{β}}}|{\mathbf{θ}} \right)$$
and

$$S\left( t|{\mathbf{θ}},{\mathbf{β}},\mathbf{x} \right) = S_{0}\left( te^{- \mathbf{x}{\mathbf{β}}}|{\mathbf{θ}} \right).$$

### Proportional Hazards Models

Proportional hazards (PH) models are defined as

$$h\left( t|{\mathbf{θ}},{\mathbf{β}},\mathbf{x} \right) = h_{0}\left( t|{\mathbf{θ}} \right)\exp\{\mathbf{x}{\mathbf{β}}\},$$
where $h_{0}\left( t|{\mathbf{θ}} \right)$ is a baseline hazard function
so that

$$f\left( t|{\mathbf{θ}},{\mathbf{β}},\mathbf{x} \right) = h_{0}\left( t|{\mathbf{θ}} \right)\exp\left\{ \mathbf{x}{\mathbf{β}} - H_{0}\left( t|{\mathbf{θ}} \right)e^{\mathbf{x}{\mathbf{β}}} \right\},$$
and

$$S\left( t|{\mathbf{θ}},{\mathbf{β}},\mathbf{x} \right) = \exp\left\{ - H_{0}\left( t|{\mathbf{θ}} \right)e^{\mathbf{x}{\mathbf{β}}} \right\}.$$

### Proportional Odds Models

Proportional Odds (PO) models are defined as

$$R\left( t|{\mathbf{θ}},{\mathbf{β}},\mathbf{x} \right) = R_{0}\left( t|{\mathbf{θ}} \right)\exp\{\mathbf{x}{\mathbf{β}}\},$$
where
$R_{0}\left( t|{\mathbf{θ}} \right) = \frac{1 - S_{0}\left( t|{\mathbf{θ}} \right)}{S_{0}\left( t|{\mathbf{θ}} \right)} = \exp\{ H_{0}\left( t|{\mathbf{θ}} \right)\} - 1$
is a baseline odds function so that

$$f\left( t|{\mathbf{θ}},{\mathbf{β}},\mathbf{x} \right) = \frac{h_{0}\left( t|{\mathbf{θ}} \right)\exp\{\mathbf{x}{\mathbf{β}} + H_{0}\left( t|{\mathbf{θ}} \right)\}}{\left\lbrack 1 + R_{0}\left( t|{\mathbf{θ}} \right)e^{\mathbf{x}{\mathbf{β}}} \right\rbrack^{2}}.$$

and

$$S\left( t|{\mathbf{θ}},{\mathbf{β}},\mathbf{x} \right) = \frac{1}{1 + R_{0}\left( t|{\mathbf{θ}} \right)e^{\mathbf{x}{\mathbf{β}}}}.$$

### Accelerated Hazard Models

Accelerated hazard (AH) models can be defined as

$$h\left( t|{\mathbf{θ}},{\mathbf{β}},\mathbf{x} \right) = h_{0}\left( t/e^{\mathbf{x}{\mathbf{β}}}|{\mathbf{θ}} \right)$$

so that

$$S\left( t|{\mathbf{θ}},{\mathbf{β}},\mathbf{x} \right) = \exp\left\{ - H_{0}\left( t/e^{\mathbf{x}{\mathbf{β}}}|{\mathbf{θ}} \right)e^{\mathbf{x}{\mathbf{β}}} \right\}$$
and
$$f\left( t|{\mathbf{θ}},{\mathbf{β}},\mathbf{x} \right) = h_{0}\left( t/e^{\mathbf{x}{\mathbf{β}}}|{\mathbf{θ}} \right)\exp\left\{ - H_{0}\left( t/e^{\mathbf{x}{\mathbf{β}}}|{\mathbf{θ}} \right)e^{\mathbf{x}{\mathbf{β}}} \right\}.$$

### Extended hazard Models

The survival function of the extended hazard (EH) model is given by:

$$S\left( t|{\mathbf{θ}},{\mathbf{β}},{\mathbf{ϕ}} \right) = \exp\left\{ - H_{0}\left( t/e^{\mathbf{x}{\mathbf{β}}}|{\mathbf{θ}} \right)\exp\left( \mathbf{x}({\mathbf{β}} + {\mathbf{ϕ}}) \right) \right\}.$$

The hazard and the probability density functions are then expressed as:

$$h\left( t|{\mathbf{θ}},{\mathbf{β}},{\mathbf{ϕ}} \right) = h_{0}\left( t/e^{\mathbf{x}{\mathbf{β}}}|{\mathbf{θ}} \right)\exp\{\mathbf{x}{\mathbf{ϕ}}\}$$
and

$$f\left( t|{\mathbf{θ}},{\mathbf{β}},{\mathbf{ϕ}} \right) = h_{0}\left( t/e^{\mathbf{x}{\mathbf{β}}}|{\mathbf{θ}} \right)\exp\{\mathbf{x}{\mathbf{β}}\}\exp\left\{ - H_{0}\left( t/e^{\mathbf{x}{\mathbf{β}}}|{\mathbf{θ}} \right)\exp\left( \mathbf{x}({\mathbf{β}} + {\mathbf{ϕ}}) \right) \right\},$$

respectively.

The EH model includes the AH, AFT and PH models as particular cases when
${\mathbf{ϕ}} = \mathbf{0}$, ${\mathbf{ϕ}} = - {\mathbf{β}}$, and
${\mathbf{β}} = \mathbf{0}$, respectively.

### Yang and Prentice Models

The survival function of the Yang and Prentice (YP) model is given by:

$$S\left( t|{\mathbf{θ}},{\mathbf{β}},{\mathbf{ϕ}} \right) = \left\lbrack 1 + \frac{\kappa_{S}}{\kappa_{L}}R_{0}\left( t|{\mathbf{θ}} \right) \right\rbrack^{- \kappa_{L}}.$$

The hazard and the probability density functions are then expressed as:

$$h\left( t|{\mathbf{θ}},{\mathbf{β}},{\mathbf{ϕ}} \right) = \frac{\kappa_{S}h_{0}\left( t|{\mathbf{θ}} \right)\exp\{ H_{0}\left( t|{\mathbf{θ}} \right)\}}{\left\lbrack 1 + \frac{\kappa_{S}}{\kappa_{L}}R_{0}\left( t|{\mathbf{θ}} \right) \right\rbrack}$$
and

$$f\left( t|{\mathbf{θ}},{\mathbf{β}},{\mathbf{ϕ}} \right) = \kappa_{S}h_{0}\left( t|{\mathbf{θ}} \right)\exp\{ H_{0}\left( t|{\mathbf{θ}} \right)\}\left\lbrack 1 + \frac{\kappa_{S}}{\kappa_{L}}R_{0}\left( t|{\mathbf{θ}} \right) \right\rbrack^{- {(1 + \kappa_{L})}},$$

respectively, where $\kappa_{S} = \exp\{\mathbf{x}{\mathbf{β}}\}$ and
$\kappa_{L} = \exp\{\mathbf{x}{\mathbf{ϕ}}\}$.

The YO model includes the PH and PO models as particular cases when
${\mathbf{ϕ}} = {\mathbf{β}}$ and ${\mathbf{ϕ}} = \mathbf{0}$,
respectively.
