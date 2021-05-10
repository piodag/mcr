\name{mc.analytical.ci}
\alias{mc.analytical.ci}
\title{Analytical Confidence Interval}
\usage{
mc.analytical.ci(b0, b1, se.b0, se.b1, n, alpha)
}
\arguments{
  \item{b0}{point estimate of intercept.}

  \item{b1}{point estimate of slope.}

  \item{se.b0}{standard error of intercept.}

  \item{se.b1}{standard error of slope.}

  \item{n}{number of observations.}

  \item{alpha}{numeric value specifying the 100(1-alpha)\%
  confidence level for the confidence interval (Default is
  0.05).}
}
\value{
2x4 matrix of estimates and confidence intervals for
intercept and slope.
}
\description{
Calculate wald confidence intervals for intercept and slope
given point estimates and standard errors.
}

