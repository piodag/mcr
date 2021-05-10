\name{mc.calc.Student}
\alias{mc.calc.Student}
\title{Student Method for Calculation of Resampling Confidence Intervals}
\usage{
mc.calc.Student(Xboot, xhat, alpha, npoints)
}
\arguments{
  \item{Xboot}{vector of point estimates for each bootstrap
  sample.  The i-th element contains the point estimate of
  the i-th bootstrap sample.}

  \item{alpha}{numeric value specifying the 100(1-alpha)\%
  confidence level for the confidence interval (Default is
  0.05).}

  \item{xhat}{global point estimate for which the
  confidence interval shall be computed.}

  \item{npoints}{number of points used for the regression
  analysis.}
}
\value{
a list with elements \item{est}{the point estimate xhat}
\item{se}{standard deviation computed from bootstrap point
estimates Xboot} \item{CI}{Confidence interval for point
estimate xhat, calculated as \eqn{xhat +/-
qt(1-alpha,n-2)*sd}.}
}
\description{
Calculate bootstrap confidence intervals for intercept,
slope or bias given a vector of bootstrap point estimates.
}
\references{
Carpenter, J., Bithell, J. (2000) Bootstrap confidence
intervals: when, which, what? A practical guide for medical
statisticians.  \emph{Stat Med}, \bold{19 (9)}, 1141--1164.
}

