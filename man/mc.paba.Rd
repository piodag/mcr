\name{mc.paba}
\alias{mc.paba}
\title{Passing-Bablok Regression}
\usage{
mc.paba(angM, X, Y, alpha = 0.05, posCor = TRUE, calcCI = TRUE,
  slope.measure = c("radian", "tangent"))
}
\arguments{
  \item{angM}{upper triangular matrix of slopes for all
  point combinations. Slopes in radian.}

  \item{X}{measurement values of reference method}

  \item{Y}{measurement values of test method}

  \item{alpha}{numeric value specifying the 100(1-alpha)%
  confidence level}

  \item{posCor}{should algorithm assume positive
  correlation, i.e. symmetry around slope 1?}

  \item{calcCI}{should confidence intervals be computed?}

  \item{slope.measure}{angular measure of pairwise slopes
  (see \code{\link{mcreg}} for details).\cr \code{"radian"} -
  for data sets with even sample numbers median slope is
  calculated as average of two central slope angles.\cr
  \code{"tangent"} - for data sets with even sample numbers
  median slope is calculated as average of two central
  slopes (tan(angle)).\cr}
}
\value{
Matrix of estimates and confidence intervals for intercept
and slope. No standard errors provided by this algorithm.
}
\description{
Passing-Bablok Regression
}

