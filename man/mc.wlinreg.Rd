\name{mc.wlinreg}
\alias{mc.wlinreg}
\title{Calculate Weighted Ordinary Linear Regression
and Estimate Standard Errors}
\usage{
mc.wlinreg(X, Y)
}
\arguments{
  \item{X}{measurement values of reference method.}

  \item{Y}{measurement values of test method.}
}
\value{
a list with elements.  \item{b0}{intercept.}
\item{b1}{slope.} \item{se.b0}{respective standard error of
intercept.} \item{se.b1}{respective standard error of
slope.} \item{xw}{weighted average of reference method
values.}
}
\description{
The weights of regression are taken as reverse squared
values of the reference method, that's why it is impossible
to achieve the calculations for zero values.
}
\references{
Neter J., Wassermann W., Kunter M.  Applied Statistical
Models.  Richard D. Irwing, INC., 1985.
}

