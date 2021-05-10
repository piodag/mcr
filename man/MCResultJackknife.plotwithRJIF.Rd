\name{MCResultJackknife.plotwithRJIF}
\alias{MCResultJackknife.plotwithRJIF}
\alias{plotwithRJIF}
\title{Plotting the Relative Jackknife Influence Function}
\usage{
MCResultJackknife.plotwithRJIF(.Object)
}
\arguments{
  \item{.Object}{object of class "MCResultJackknife" or
  "MCResultResampling"}
}
\description{
The function draws reference method vs. test method as
scatter plot. Observations with high influence (relative
jackknife influence function is greater than 2) are
highlighted as red points.
}
\examples{
#library("mcr")
    data(creatinine,package="mcr")
    x <- creatinine$serum.crea
    y <- creatinine$plasma.crea
    # Deming regression fit.
    # The confidence intervals for regression coefficients
    # are calculated with jackknife method
    model <- mcreg( x,y,error.ratio=1,method.reg="Deming", method.ci="jackknife",
                     mref.name = "serum.crea", mtest.name = "plasma.crea", na.rm=TRUE )
    plotwithRJIF(model)
}
\references{
Efron, B. (1990) Jackknife-After-Bootstrap Standard Errors
and Influence Functions.  Technical Report , \bold{N 134}.
}

