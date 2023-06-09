% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/residuals.R
\name{ggresiduals.survstan}
\alias{ggresiduals.survstan}
\title{ggresiduals method for survstan models}
\usage{
\method{ggresiduals}{survstan}(object, type = c("coxsnell", "martingale", "deviance"), ...)
}
\arguments{
\item{object}{a fitted model object of the class survstan.}

\item{type}{type of residuals used in the plot: coxsnell (default), martingale and deviance.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
the desired residual plot.
}
\description{
ggresiduals method for survstan models
}
\details{
This function produces residuals plots of Cox-Snell residuals, martingale residuals and deviance residuals.
}
\examples{
\donttest{
library(survstan)
ovarian$rx <- as.factor(ovarian$rx)
fit <- aftreg(Surv(futime, fustat) ~ age + rx, data = ovarian, baseline = "weibull", init = 0)
ggresiduals(fit, type = "coxsnell")
ggresiduals(fit, type = "martingale")
ggresiduals(fit, type = "deviance")
}

}
