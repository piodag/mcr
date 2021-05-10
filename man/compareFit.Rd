\name{compareFit}
\alias{compareFit}
\title{Graphical Comparison of Regression Parameters and Associated Confidence Intervals}
\usage{
compareFit(...)
}
\arguments{
  \item{...}{list of fitted models, i.e. objects of
  "MCResult" sub-classes.}
}
\description{
Graphical comparison of regression parameters (intercept
and slope) and their associated 100(1-alpha)\% confidence
intervals for multiple fitted models of 'MCResult'
sub-classes.
}
\examples{
library("mcr")
     data("creatinine", package="mcr")
     fit.lr <- mcreg(as.matrix(creatinine), method.reg="LinReg", na.rm=TRUE)
     fit.wlr <- mcreg(as.matrix(creatinine), method.reg="WLinReg", na.rm=TRUE)
     compareFit( fit.lr, fit.wlr )
}

