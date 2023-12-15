# survstan 0.0.1

* Added a `NEWS.md` file to track changes to the package.

# survstan 0.0.2

- Implementation of the Yang and Prentice (YP) model for dealing with survival data with possibly crossing survival curves.
- Implementation of survfit.survstan method for all implemented models.
- Improvement in the output of the anova.survstan function: now the chosen regression model and baseline distribution are displayed in the ANOVA table.
- Addition of dist argument in the model fitting functions (aftreg(), ahreg(), phreg(), poreg() and ypreg()) as an alternative way to specify the baseline distribution (for compability with the survival::survreg() function).
- Addition of init argument in the model fitting functions (aftreg(), ahreg(), phreg(), poreg() and ypreg()) with default value equals to 0 to fix the initial starting values when calling the rstan::optimizing() function.
- Addition of extractAIC.survstan method required in the stats::step() function for carrying out model selection via stepwise procedures.
- Addition of logLik.survstan method to extract the log-likelihood value of a fitted model.
- AIC.survstan method now returns a data.frame with the aic values along with the corresponding number of model parameters when more than one fit is passed to the method.
- Addition of gastric and ipass datasets.

# survstan 0.0.3

- The Weibull distribution is now the default baseline distribution (replacing the exponential distribution in previous versions).
- Implementation of Birnbaum-Saunders (fatigue) distribution.
- Implementation of the function rank_models to rank survstan models according to their AICs.


# survstan 0.0.4

- The survstan package now requires rstan version 2.26.


# survstan 0.0.5

- Inclusion of emmeans method for survstan models.
- Implementation of gamma and Rayleigh distributions.
- Implementation of se method for computation of standard errors.
- Correction of covariance matrix for positive parameters.
- Application of delta method to avoid negative lower bounds in confidence intervals for positive parameters.
- Bug correction in the implementation of survival functions.
