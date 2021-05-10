\name{MCResult.getResiduals}
\alias{MCResult.getResiduals}
\alias{getResiduals}
\title{Get Regression Residuals}
\usage{
MCResult.getResiduals(.Object)
}
\arguments{
  \item{.Object}{object of class "MCResult".}
}
\value{
residuals as data frame.
}
\description{
This function returnes residuals in x-direction (x-xhat),
in y-direction(y-yhat) and optimized residuals. The
optimized residuals correspond to distances between data
points and the regression line which were optimized for
regression coefficients estimation. In case of
Passing-Bablok Regression orthogonal residuals will be
returned as optimized residuals . The residuals in
x-direction are interesting for regression types which
assume errors in both variables (deming, weighted deming,
Passing-Bablok), particularily for checking of model
assumptions.
}
\seealso{
\code{\link{plotResiduals}}
}

