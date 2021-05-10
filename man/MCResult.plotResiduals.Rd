\name{MCResult.plotResiduals}
\alias{MCResult.plotResiduals}
\alias{plotResiduals}
\title{Plot Residuals of an MCResult Object}
\usage{
MCResult.plotResiduals(.Object, res.type = c("optimized", "y", "x"),
  xaxis = c("yhat", "both", "xhat"), ref.line = TRUE,
  ref.line.col = "red", ref.line.lty = 2, ref.line.lwd = 1, main = NULL,
  xlab = NULL, ylab = NULL, add.grid = TRUE, ...)
}
\arguments{
  \item{.Object}{object of type "MCResult".}

  \item{res.type}{If res.type="y" the difference between
  the test method and it's prediction will be drawn. If
  res.type="x" the reference method and it's prediction
  will be drawn. In case ordinary and weighted ordinary
  linear regression this difference will be zero.}

  \item{xaxis}{Values on the x-axis. One can choose from
  estimated values of x (xaxis=\code{"xhat"}), y
  (xaxis=\code{"xhat"}) or the mean of estimated values of
  x and y (\code{xaxis="both"}). If res.type="optimized"
  the proper type of residuals for each regression will be
  drawn.}

  \item{ref.line}{logical value. If \code{ref.line = TRUE}
  (default), the reference line will be drawn.}

  \item{ref.line.col}{reference line color.}

  \item{ref.line.lty}{reference line type.}

  \item{ref.line.lwd}{reference line width.}

  \item{xlab}{label for the x-axis}

  \item{ylab}{label for the y-axis}

  \item{add.grid}{logical value. If \code{add.grid = TRUE}
  (default) the gridlines will be drawn.}

  \item{main}{character string specifying the main title of
  the plot}

  \item{...}{further graphical parameters}
}
\description{
Plot Residuals of an MCResult Object
}
\examples{
data(creatinine,package="mcr")
    x <- creatinine$serum.crea
    y <- creatinine$plasma.crea

    # Deming regression fit.
    # The confidence intercals for regression coefficients
    # are calculated with analytical method
    model <- mcreg( x,y,error.ratio=1,method.reg="WDeming", method.ci="jackknife",
                     mref.name = "serum.crea", mtest.name = "plasma.crea", na.rm=TRUE )
    plotResiduals(model, res.type="optimized", xaxis="both" )
    plotResiduals(model, res.type="y", xaxis="yhat")
}
\seealso{
\code{\link{getResiduals}}, \code{\link{plot.mcr}},
\code{\link{plotDifference}}, \code{\link{plotBias}},
\code{\link{compareFit}}
}

