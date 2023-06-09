% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aic.R
\name{AIC.survstan}
\alias{AIC.survstan}
\title{Akaike information criterion}
\usage{
\method{AIC}{survstan}(object, ..., k = 2)
}
\arguments{
\item{object}{an object of the class survstan.}

\item{...}{further arguments passed to or from other methods.}

\item{k}{numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.}
}
\value{
the Akaike information criterion value when a single model is passed to the function; otherwise, a data.frame with the Akaike information criterion values and the number of parameters is returned.
}
\description{
Akaike information criterion
}
\examples{
\donttest{
library(survstan)
fit1 <- aftreg(Surv(futime, fustat) ~ 1, data = ovarian, baseline = "weibull", init = 0)
fit2 <- aftreg(Surv(futime, fustat) ~ rx, data = ovarian, baseline = "weibull", init = 0)
fit3 <- aftreg(Surv(futime, fustat) ~ ecog.ps + rx, data = ovarian, baseline = "weibull", init = 0)
AIC(fit1, fit2, fit3)
}

}
