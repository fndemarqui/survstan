% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/piecewise.R
\name{time_grid}
\alias{time_grid}
\title{Time grid}
\usage{
time_grid(time, event, m = NULL)
}
\arguments{
\item{time}{Vector of failure times}

\item{event}{Vector of failure indicators}

\item{m}{Optional. Number of intervals. If \code{NULL}, the number of intervals is
set to be equal to the number of distinct observed failure times.}
}
\value{
Time grid.
}
\description{
Time grid
}
