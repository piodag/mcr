\name{mc.deming}
\alias{mc.deming}
\title{Calculate Unweighted Deming Regression and Estimate Standard Errors}
\usage{
mc.deming(X, Y, error.ratio)
}
\arguments{
  \item{X}{measurement values of reference method.}

  \item{Y}{measurement values of test method.}

  \item{error.ratio}{ratio of measurement error of
  reference method to measurement error of test method.}
}
\value{
a list with elements \item{b0}{intercept.}
\item{b1}{slope.} \item{se.b0}{respective standard error of
intercept.} \item{se.b1}{respective standard error of
slope.} \item{xw}{average of reference method values.}
}
\description{
Calculate Unweighted Deming Regression and Estimate
Standard Errors
}
\references{
Linnet K.  Evaluation of Regression Procedures for Methods
Comparison Studies.  CLIN. CHEM. 39/3, 424-432 (1993).

Linnet K.  Estimation of the Linear Relationship between
the Measurements of two Methods with Proportional Errors.
STATISTICS IN MEDICINE, Vol. 9, 1463-1473 (1990).
}

