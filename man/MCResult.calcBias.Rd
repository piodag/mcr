\name{MCResult.calcBias}
\alias{MCResult.calcBias}
\alias{calcBias}
\title{Systematical Bias Between Reference Method and Test Method}
\usage{
MCResult.calcBias(.Object, x.levels, type = c("absolute", "proportional"),
  percent = TRUE, alpha = 0.05, ...)
}
\arguments{
  \item{.Object}{object of class "MCResult".}

  \item{type}{One can choose between absolute (default) and
  proportional bias (\code{Bias(Xc)/Xc}).}

  \item{percent}{logical value. If \code{percent = TRUE}
  the proportional bias will be calculated in percent.}

  \item{x.levels}{a numeric vector with decision points for
  which bias schould be calculated.}

  \item{alpha}{numeric value specifying the
  100(1-\code{alpha})\% confidence level of the confidence
  interval (Default is 0.05).}

  \item{...}{further parameters}
}
\value{
response and corresponding confidence interval for each
decision point from x.levels.
}
\description{
Calculate systematical bias between reference and test
methods at the decision point Xc as \eqn{ Bias(Xc) =
Intercept + (Slope-1) * Xc} with corresponding confidence
intervals.
}
\examples{
#library("mcr")
    data(creatinine,package="mcr")
    x <- creatinine$serum.crea
    y <- creatinine$plasma.crea

    # Deming regression fit.
    # The confidence intervals for regression coefficients
    # are calculated with analytical method
    model <- mcreg( x,y,error.ratio = 1,method.reg = "Deming", method.ci = "analytical",

                     mref.name = "serum.crea", mtest.name = "plasma.crea", na.rm=TRUE )
    # Now we calculate the systematical bias
    # between the testmethod and the reference method
    # at the medical decision points 1, 2 and 3

    calcBias( model, x.levels = c(1,2,3))
    calcBias( model, x.levels = c(1,2,3), type = "proportional")
    calcBias( model, x.levels = c(1,2,3), type = "proportional", percent = FALSE)
}
\seealso{
\code{\link{plotBias}}
}

