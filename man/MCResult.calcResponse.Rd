\name{MCResult.calcResponse}
\alias{MCResult.calcResponse}
\alias{calcResponse}
\title{Calculate Response with Confidence Interval.}
\usage{
MCResult.calcResponse(.Object, x.levels, alpha, ...)
}
\arguments{
  \item{.Object}{object of class "MCResult".}

  \item{x.levels}{a numeric vector with points for which
  response schould be calculated.}

  \item{alpha}{numeric value specifying the
  100(1-\code{alpha})\% confidence level of the confidence
  interval (Default is 0.05).}

  \item{...}{further parameters}
}
\value{
response and corresponding confidence interval for each
point in vector \code{x.levels}.
}
\description{
Calculate Response \eqn{Intercept + Slope * Refrencemethod}
with Corresponding Confidence Interval
}
\examples{
#library("mcr")
    data(creatinine,package="mcr")
    x <- creatinine$serum.crea
    y <- creatinine$plasma.crea
    # Deming regression fit.
    # The confidence intercals for regression coefficients
    # are calculated with analytical method
    model <- mcreg( x,y,error.ratio=1,method.reg="Deming", method.ci="analytical",
                     mref.name = "serum.crea", mtest.name = "plasma.crea", na.rm=TRUE )
    calcResponse(model, x.levels=c(1,2,3))
}
\seealso{
\code{\link{calcBias}}
}

