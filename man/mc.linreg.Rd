\name{mc.linreg}
\alias{mc.linreg}
\title{Calculate ordinary linear Regression
and Estimate Standard Errors}
\usage{
mc.linreg(X, Y)
}
\arguments{
  \item{X}{measurement values of reference method.}

  \item{Y}{measurement values of test method.}
}
\value{
a list with elements \item{b0}{intercept.}
\item{b1}{slope.} \item{se.b0}{respective standard error of
intercept.} \item{se.b1}{respective standard error of
slope.} \item{xw}{average of reference method values.}
}
\description{
Calculate ordinary linear Regression and Estimate Standard
Errors
}
\references{
Neter J., Wassermann W., Kunter M.  Applied Statistical
Models.  Richard D. Irwing, INC., 1985.
}

