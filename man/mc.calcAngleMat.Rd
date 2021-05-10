\name{mc.calcAngleMat}
\alias{mc.calcAngleMat}
\title{Calculate Matrix of All Pair-wise Slope Angles}
\usage{
mc.calcAngleMat(X, Y, posCor = TRUE)
}
\arguments{
  \item{X}{measurement values of reference method.}

  \item{Y}{measurement values of test method.}

  \item{posCor}{should algorithm assume positive
  correlation, i.e. symmetry around slope 1?}
}
\value{
Upper triangular matrix of slopes for all point
combinations. Slopes in radian.
}
\description{
This version is implemented in C for computational
efficiency.
}

