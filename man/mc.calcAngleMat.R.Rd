\name{mc.calcAngleMat.R}
\alias{mc.calcAngleMat.R}
\title{Calculate Matrix of All Pair-wise Slope Angles}
\usage{
mc.calcAngleMat.R(X, Y, posCor = TRUE)
}
\arguments{
  \item{X}{measurement values of reference method.}

  \item{Y}{measurement values of test method.}

  \item{posCor}{should the algorithm assume positive
  correlation, i.e. symmetry around slope 1?}
}
\value{
Upper triangular matrix of slopes for all point
combinations. Slopes in radian.
}
\description{
This is a very slow R version. It should not be called
except for debugging purposes.
}

