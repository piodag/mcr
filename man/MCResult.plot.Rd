\name{MCResult.plot}
\alias{MCResult.plot}
\alias{plot.mcr}
\title{Scatter Plot Method X vs. Method Y}
\usage{
MCResult.plot(x, alpha = 0.05, xn = 20, equal.axis = FALSE, xlim = NULL,
  ylim = NULL, x.lab = x@mnames[1], y.lab = x@mnames[2], add = FALSE,
  draw.points = TRUE, points.col = "black", points.pch = 1,
  points.cex = 0.8, reg = TRUE, reg.col = NULL, reg.lty = 1,
  reg.lwd = 2, identity = TRUE, identity.col = NULL, identity.lty = 2,
  identity.lwd = 1, ci.area = TRUE, ci.area.col = NULL,
  ci.border = FALSE, ci.border.col = NULL, ci.border.lty = 2,
  ci.border.lwd = 1, add.legend = TRUE, legend.place = c("topleft",
  "topright", "bottomleft", "bottomright"), main = NULL, sub = NULL,
  add.cor = TRUE, cor.method = c("pearson", "kendall", "spearman"),
  add.grid = TRUE, ...)
}
\arguments{
  \item{x}{object of class "MCResult".}

  \item{alpha}{numeric value specifying the
  100(1-\code{alpha})\% confidence bounds.}

  \item{xn}{number of points (default 20) for calculation
  of confidence bounds.}

  \item{draw.points}{logical value. If
  \code{draw.points=TRUE}, the data points will be drawn.}

  \item{xlim}{limits of the x-axis. If \code{xlim=NULL} the
  x-limits will be calculated automatically.}

  \item{ylim}{limits of the y-axis. If \code{ylim=NULL} the
  y-limits will be calculated automatically.}

  \item{x.lab}{label of x-axis. Default is the name of
  reference method.}

  \item{y.lab}{label of y-axis. Default is the name of test
  method.}

  \item{equal.axis}{logical value. If
  \code{equal.axis=TRUE} x-axis will be equal to y-axis.}

  \item{add}{logical value. If \code{add=TRUE}, the plot
  will be drawn in current graphical window.}

  \item{points.col}{Color of data points.}

  \item{points.pch}{Type of data points (see
  \code{par()}).}

  \item{points.cex}{Size of data points (see
  \code{par()}).}

  \item{reg}{Logical value. If \code{reg=TRUE}, the
  regression line will be drawn.}

  \item{reg.col}{Color of regression line.}

  \item{reg.lty}{Type of regression line.}

  \item{reg.lwd}{The width of regression line.}

  \item{identity}{logical value. If \code{identity=TRUE}
  the identity line will be drawn.}

  \item{identity.col}{The color of identity line.}

  \item{identity.lty}{The type of identity line.}

  \item{identity.lwd}{the width of identity line.}

  \item{ci.area}{logical value. If \code{ci.area=TRUE}
  (default) the confidence area will be drawn.}

  \item{ci.area.col}{the color of confidence area.}

  \item{ci.border}{logical value. If \code{ci.border=TRUE}
  the confidence limits will be drawn.}

  \item{ci.border.col}{The color of confidence limits.}

  \item{ci.border.lty}{The line type of confidence limits.}

  \item{ci.border.lwd}{The line width of confidence
  limits.}

  \item{add.legend}{logical value. If
  \code{add.legend=FALSE} the plot will not have any
  legend.}

  \item{legend.place}{The position of legend:
  "topleft","topright","bottomleft","bottomright".}

  \item{main}{String value. The main title of plot. If
  \code{main=NULL} it will include regression name.}

  \item{sub}{String value. The subtitle of plot. If
  \code{sub=NULL} and \code{ci.border=TRUE} or
  \code{ci.area=TRUE} it will include the art of confidence
  bounds calculation.}

  \item{add.cor}{Logical value. If \code{add.cor=TRUE} the
  correlation coefficient will be shown.}

  \item{cor.method}{a character string indicating which
  correlation coefficient is to be computed. One of
  "pearson" (default), "kendall", or "spearman", can be
  abbreviated.}

  \item{add.grid}{Logical value. If \code{add.grid=TRUE}
  (default) the gridlines will be drawn.}

  \item{...}{further graphical parameters}
}
\description{
Plot method X (reference) vs. method Y (test) with
(optional) line of identity, regression line and confidence
bounds for response.
}
\examples{
library(mcr)
 data(creatinine,package="mcr")
 creatinine <- creatinine[complete.cases(creatinine),]
  x <- creatinine$serum.crea
  y <- creatinine$plasma.crea

  m1 <- mcreg(x,y,method.reg="Deming",  mref.name="serum.crea",
                                        mtest.name="plasma.crea", na.rm=TRUE)
  m2 <- mcreg(x,y,method.reg="WDeming", method.ci="jackknife",
                                        mref.name="serum.crea",
                                        mtest.name="plasma.crea", na.rm=TRUE)

  plot(m1,  xlim=c(0.5,3),ylim=c(0.5,3), add.legend=FALSE,
                           main="Deming vs. weighted Deming regression",
                           points.pch=19,ci.area=TRUE, ci.area.col=grey(0.9),
                           identity=FALSE, add.grid=FALSE, sub="")
  plot(m2, ci.area=FALSE, ci.border=TRUE, ci.border.col="red3",
                           reg.col="red3", add.legend=FALSE,
                           draw.points=FALSE,add=TRUE)

  includeLegend(place="topleft",models=list(m1,m2),
                           colors=c("darkblue","red"), design="1", digits=2)
}
\seealso{
\code{\link{plotBias}}, \code{\link{plotResiduals}},
\code{\link{plotDifference}},
\code{\link{compareFit}},\code{\link{includeLegend}}
}

