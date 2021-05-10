\name{mc.mdemingConstCV}
\alias{mc.mdemingConstCV}
\title{Calculate M Weighted Deming Regression}
\usage{
mc.mdemingConstCV(X, Y, error.ratio, iter.max = 30, threshold = 1e-06)
}
\arguments{
  \item{X}{measurement values of reference method.}

  \item{Y}{measurement values of test method.}

  \item{error.ratio}{ratio between squared measurement
  errors of reference- and test method, necessary for
  Deming regression (Default is 1).}

  \item{iter.max}{maximal number of iterations.}

  \item{threshold}{threshold value.}
}
\value{
a list with elements \item{b0}{intercept.}
\item{b1}{slope.} \item{xw}{average of reference method
values.} \item{iter}{number of iterations.}
}
\description{
Calculate M-weighted deming regression with iterative
algorithm suggested by Linnet. There is no
guarantee that the algorithm always converges.
}
\references{
Linnet K.  Evaluation of Regression Procedures for Methods
Comparison Studies.  CLIN. CHEM. 39/3, 424-432 (1993).

Linnet K.  Estimation of the Linear Relationship between
the Measurements of two Methods with Proportional Errors.
STATISTICS IN MEDICINE, Vol. 9, 1463-1473 (1990).
}

