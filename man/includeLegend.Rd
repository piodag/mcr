\name{includeLegend}
\alias{includeLegend}
\title{Include Legend}
\usage{
includeLegend(models = list(), digits = 2, design = paste(1:2),
  place = c("topleft", "topright", "bottomleft", "bottomright"), colors,
  lty = rep(1, length(models)), lwd = rep(2, length(models)),
  box.lty = "blank", cex = 0.8, bg = "white", inset = c(0.01, 0.01),
  bias = FALSE, model.names = NULL, ...)
}
\arguments{
  \item{models}{list of length n with Objects of class
  "MCResult".}

  \item{digits}{number of digits in Coefficients.}

  \item{design}{type of legend design. There are two
  possible designs: "1" and "2" (See example).}

  \item{place}{place for Legend:
  "topleft","topright","bottomleft" or "bottomright".}

  \item{colors}{vector of length n with color of regression
  lines.}

  \item{lty}{vector of length n with type of regression
  lines.}

  \item{lwd}{vector of length n with thickness of
  regression lines.}

  \item{box.lty}{box line-type}

  \item{cex}{numeric value representing the plotting symbol
  magnification factor}

  \item{bg}{the background-color of the legend box}

  \item{inset}{inset distance(s) from the margins as a
  fraction of the plot region when legend is placed by
  keyword.}

  \item{bias}{logical value. If bias = TRUE, it will be
  drawn a legend for \code{plotBias()} function.}

  \item{model.names}{legend names for different models. If
  NULL the regression type will be used.}

  \item{...}{other parameters of function legend().}
}
\value{
Legend in plot.
}
\description{
Include legend in regression plot (function \code{plot()})
or in bias plot (function \code{plotBias ()}) with two or
more lines.
}
\examples{
#library("mcr")

 data(creatinine,package="mcr")
 x <- creatinine$serum.crea
 y <- creatinine$plasma.crea

 m1 <- mcreg(x,y,method.reg="Deming", mref.name="serum.crea",
                                        mtest.name="plasma.crea", na.rm=TRUE)
 m2 <- mcreg(x,y,method.reg="WDeming", method.ci="jackknife",
                                         mref.name="serum.crea",
                                         mtest.name="plasma.crea", na.rm=TRUE)

 plot(m1,  XLIM=c(0.5,3),YLIM=c(0.5,3), Legend=FALSE,
                          Title="Deming vs. weighted Deming regression",
                          Points.pch=19,ci.area=TRUE, ci.area.col=grey(0.9),
                          identity=FALSE, Grid=FALSE, Sub="")
 plot(m2, ci.area=FALSE, ci.border=TRUE, ci.border.col="red3",
                          reg.col="red3", Legend=FALSE,add=TRUE,
                          Points=FALSE, identity=FALSE, Grid=FALSE)

 includeLegend(place="topleft",models=list(m1,m2),
                          colors=c("darkblue","red"), design="1", digits=2)
}
\seealso{
\code{\link{plot.mcr}}, \code{\link{plotBias}},
\code{\link{plotResiduals}}, \code{\link{plotDifference}},
\code{\link{compareFit}}
}

