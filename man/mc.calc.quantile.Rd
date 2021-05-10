\name{mc.calc.quantile}
\alias{mc.calc.quantile}
\title{Quantile Method for Calculation of Resampling Confidence Intervals}
\usage{
mc.calc.quantile(Xboot, alpha)
}
\arguments{
  \item{Xboot}{vector of point estimates for bootstrap
  samples.  The i-th element contains point estimate of the
  i-th bootstrap sample.}

  \item{alpha}{numeric value specifying the 100(1-alpha)\%
  confidence level for the confidence interval (Default is
  0.05).}
}
\value{
a list with elements \item{est}{median of bootstrap point
estimates Xboot.} \item{CI}{confidence interval for point
estimate 'est', calculated as quantiles.}
}
\description{
Calculate bootstrap confidence intervals for intercept,
slope or bias given the vector of bootstrap point
estimates.
}
\references{
B. Efron and RJ. Tibshirani (1994) An Introduction to the
Bootstrap.  \emph{Chapman & Hall}.
}

