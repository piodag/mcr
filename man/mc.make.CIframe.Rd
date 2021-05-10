\name{mc.make.CIframe}
\alias{mc.make.CIframe}
\title{Returns Results of Calculations in Matrix Form}
\usage{
mc.make.CIframe(b0, b1, se.b0, se.b1, CI.b0, CI.b1)
}
\arguments{
  \item{b0}{point estimate for intercept.}

  \item{b1}{point estimate for slope.}

  \item{se.b0}{standard error of intercept estimate.}

  \item{se.b1}{standard error of slope estimate.}

  \item{CI.b0}{numeric vector of length 2 - confidence
  interval for intercept.}

  \item{CI.b1}{numeric vector of length 2 - confidence
  interval for slope.}
}
\value{
2x4 matrix of estimates and confidence intervals for
intercept and slope.
}
\description{
Returns Results of Calculations in Matrix Form
}

