\name{mc.calcTstar}
\alias{mc.calcTstar}
\title{Compute Resampling T-statistic.}
\usage{
mc.calcTstar(.Object, x.levels, iter.max = 30, threshold = 1e-06)
}
\arguments{
  \item{.Object}{object of class "MCResultResampling".}

  \item{x.levels}{a numeic vector of clinical desision
  points of interest.}

  \item{iter.max}{maximal number of iterations for
  calculation of weighted deming regression.}

  \item{threshold}{threshold for calculation of weighted
  deming regression.}
}
\value{
Tstar numeric vector containing resampling pivot statistic.
}
\description{
Compute Resampling T-statistic. for Calculation of
t-Bootstrap Confidence Intervals.
}
\references{
Carpenter J., Bithell J.  Bootstrap confidence intervals:
when, which, what? A practical guide for medical
statisticians.  Stat Med, 19 (9), 1141-1164 (2000).
}

