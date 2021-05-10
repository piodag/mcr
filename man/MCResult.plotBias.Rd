\name{MCResult.plotBias}
\alias{MCResult.plotBias}
\alias{plotBias}
\title{Plot Estimated Systematical Bias with Confidence Bounds}
\usage{
MCResult.plotBias(x, xn = 100, alpha = 0.05, add = FALSE, prop = FALSE,
  xlim = NULL, ylim = NULL, bias = TRUE, bias.lty = 1, bias.lwd = 2,
  bias.col = NULL, ci.area = TRUE, ci.area.col = NULL,
  ci.border = FALSE, ci.border.col = NULL, ci.border.lwd = 1,
  ci.border.lty = 2, zeroline = TRUE, zeroline.col = NULL,
  zeroline.lty = 2, zeroline.lwd = 1, main = NULL, sub = NULL,
  add.grid = TRUE, xlab = NULL, ylab = NULL, cut.point = NULL,
  cut.point.col = "red", cut.point.lwd = 2, cut.point.lty = 1, ...)
}
\arguments{
  \item{x}{object of class "MCResult".}

  \item{xn}{# number of poits for drawing of confidence
  bounds/area.}

  \item{add}{logical value. If \code{add=TRUE}, the grafic
  will be drawn in current grafical window.}

  \item{prop}{a logical value. If \code{prop=TRUE} the
  proportional bias \eqn{ % bias(Xc) = [ Intercept +
  (Slope-1) * Xc ] / Xc} will be drawn.}

  \item{xlim}{limits of the x-axis. If \code{xlim=NULL} the
  x-limits will be calculated automatically.}

  \item{ylim}{limits of the y-axis. If \code{ylim=NULL} the
  y-limits will be calculated automatically.}

  \item{bias}{logical value. If \code{identity=TRUE} the
  bias line will be drawn. If ci.bounds=FALSE and
  ci.area=FALSE the bias line will be drawn always.}

  \item{bias.col}{color of the bias line.}

  \item{bias.lty}{type of the bias line.}

  \item{bias.lwd}{width of the bias line.}

  \item{zeroline}{logical value. If \code{zeroline=TRUE}
  the zero-line will be drawn.}

  \item{zeroline.col}{color of the zero-line.}

  \item{zeroline.lty}{type of the zero-line.}

  \item{zeroline.lwd}{width of the zero-line.}

  \item{ci.area}{logical value. If \code{ci.area=TRUE}
  (default) the confidence area will be drawn.}

  \item{ci.border}{logical value. If ci.border=TRUE the
  confidence limits will be drawn.}

  \item{ci.area.col}{color of the confidence area.}

  \item{ci.border.col}{color of the confidence limits.}

  \item{ci.border.lty}{line type of confidence limits.}

  \item{ci.border.lwd}{line width of confidence limits.}

  \item{cut.point}{numeric value. Decision level of
  interest.}

  \item{cut.point.col}{color of the confidence bounds at
  the required decision level.}

  \item{cut.point.lty}{line type of the confidence bounds
  at the required decision level.}

  \item{cut.point.lwd}{line width of the confidence bounds
  at the required decision level.}

  \item{main}{character string. The main title of plot. If
  \code{main = NULL} it will include regression name.}

  \item{sub}{character string. The subtitle of plot. If
  \code{sub=NULL} and \code{ci.border=TRUE} or
  \code{ci.area=TRUE} it will include the art of confidence
  bounds calculation.}

  \item{add.grid}{logical value. If \code{grid=TRUE}
  (default) the gridlines will be drawn.}

  \item{xlab}{label for the x-axis}

  \item{ylab}{label for the y-axis}

  \item{alpha}{numeric value specifying the
  100(1-\code{alpha})\% confidence level of confidence
  intervals (Default is 0.05).}

  \item{...}{further graphical parameters}
}
\description{
This function plots the estimated systematical bias \eqn{(
Intercept + Slope * Refrencemethod ) - Referencemethod}
with confidence bounds, covering the whole range of
reference method X or only part of it.
}
\examples{
#library("mcr")
data(creatinine,package="mcr")

creatinine <- creatinine[complete.cases(creatinine),]
x <- creatinine$serum.crea
y <- creatinine$plasma.crea

# Calculation of models
m1 <- mcreg(x,y,method.reg="WDeming", method.ci="jackknife",
                mref.name="serum.crea",mtest.name="plasma.crea", na.rm=TRUE)
m2 <- mcreg(x,y,method.reg="WDeming", method.ci="bootstrap",
                method.bootstrap.ci="BCa",mref.name="serum.crea",
                mtest.name="plasma.crea", na.rm=TRUE)

# Grafical comparison of systematical Bias of two models
plotBias(m1, zeroline=TRUE,zeroline.col="black",zeroline.lty=1,
                ci.area=TRUE,ci.border=FALSE, ci.area.col=grey(0.9),
                main = "Bias between serum and plasma creatinine",
                sub="Comparison of Jackknife and BCa-Bootstrap confidence bounds ")
plotBias(m2, ci.area=FALSE, ci.border=TRUE, ci.border.lwd=2,
                ci.border.col="red",bias=FALSE ,add=TRUE)
includeLegend(place="topleft",models=list(m1,m2), lwd=c(10,2),
                lty=c(2,1),colors=c(grey(0.9),"red"), bias=TRUE,
                design="1", digits=4)

# Drawing of proportional bias
plotBias(m1, ci.area=FALSE, ci.border=TRUE)
plotBias(m1, ci.area=FALSE, ci.border=TRUE, prop=TRUE)
plotBias(m1, ci.area=FALSE, ci.border=TRUE, prop=TRUE, cut.point=0.6)
plotBias(m1, ci.area=FALSE, ci.border=TRUE, prop=TRUE, cut.point=0.6,
             xlim=c(0.4,0.8),cut.point.col="orange", cut.point.lwd=3, main ="")
}
\seealso{
\code{\link{calcBias}}, \code{\link{plot.mcr}},
\code{\link{plotResiduals}}, \code{\link{plotDifference}},
\code{\link{compareFit}}
}

