\name{mc.calc.bca}
\alias{mc.calc.bca}
\title{Bias Corrected and Accelerated Resampling Confidence Interval}
\usage{
mc.calc.bca(Xboot, Xjack, xhat, alpha)
}
\arguments{
  \item{Xboot}{vector of point estimates for bootstrap
  samples.  The i-th element contains point estimate of the
  i-th bootstrap sample.}

  \item{Xjack}{vector of point estimates for jackknife
  samples.  The i-th element contains point estimate of the
  dataset without i-th observation.}

  \item{xhat}{point estimate for the complete data set
  (scalar).}

  \item{alpha}{numeric value specifying the 100(1-alpha)\%
  confidence level for the confidence interval (Default is
  0.05).}
}
\value{
a list with elements \item{est}{point estimate for the
complete data set (xhat).} \item{CI}{confidence interval
for point estimate.}
}
\description{
Calculate resampling BCa confidence intervals for
intercept, slope or bias given a vector of bootstrap and
jackknife point estimates.
}
\references{
Carpenter, J., Bithell, J. (2000) Bootstrap confidence
intervals: when, which, what? A practical guide for medical
statisticians.  \emph{Stat Med}, \bold{19 (9)}, 1141--1164.
}

