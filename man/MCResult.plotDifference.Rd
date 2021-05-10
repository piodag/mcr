\name{MCResult.plotDifference}
\alias{MCResult.plotDifference}
\alias{plotDifference}
\title{Bland-Altman Plot}
\usage{
MCResult.plotDifference(.Object, xlab = NULL, ylab = NULL,
  ref.line = TRUE, ref.line.col = "black", ref.line.lty = 1,
  ref.line.lwd = 1, bias.line.lty = 1, bias.line.lwd = 1,
  bias.line.col = "red", bias.text.col = NULL, bias.text.cex = 0.8,
  loa.line.lty = 2, loa.line.lwd = 1, loa.line.col = "red",
  loa.text.col = NULL, plot.type = 3, main = NULL, cex = 0.8,
  digits = 2, add.grid = TRUE, ylim = NULL, ...)
}
\arguments{
  \item{.Object}{object of class "MCResult".}

  \item{xlab}{label for the x-axis}

  \item{ylab}{label for the y-axis}

  \item{digits}{number of decimal places for the difference
  of means and standard deviation appearing in the plot.}

  \item{plot.type}{integer specifying a specific
  Bland-Altman plot modification (default is 3).  Possible
  choices are: 1 - difference plot X vs. Y-X with null-line
  and mean plus confidence intervals.\cr 2 - difference
  plot X vs. (Y-X)/X (relative differences) with null-line
  and mean.\cr 3 - difference plot 0.5*(X+Y) vs. Y-X with
  null-line and mean plus confidence intervals.\cr 4 -
  difference plot 0.5*(X+Y) vs. (Y-X)/X (relative
  differences) with null-line.\cr 5 - difference plot
  rank(X) vs. Y-X with null-line and mean plus confidence
  intervals.\cr 6 - difference plot rank(X) vs. (Y-X)/X
  (relative differences) with null-line and mean.\cr 7 -
  difference plot sqrt(X*Y) vs. Y/X with null-line and mean
  plus confidence intervals calculated with help of
  log-transformation.\cr 8 - difference plot 0.5*(X+Y) vs.
  (Y-X) / (0.5*(X+Y)) with null-line.\cr}

  \item{main}{plot title.}

  \item{ref.line}{logical value. If \code{ref.line=TRUE}
  (default), the reference line will be drawn.}

  \item{ref.line.col}{reference line color.}

  \item{ref.line.lty}{reference line type.}

  \item{ref.line.lwd}{reference line width.}

  \item{bias.line.lty}{line type for estimated bias.}

  \item{bias.line.lwd}{line width for estimated bias.}

  \item{bias.line.col}{color of the line for estimated
  bias.}

  \item{bias.text.col}{color of the label for estimated
  bias (defaults to the same as \code{bias.line.col}.)}

  \item{bias.text.cex}{The magnification to be used for the
  label for estimated bias}

  \item{loa.line.lty}{line type for estimated limits of
  agreement.}

  \item{loa.line.lwd}{line width for estimated limits of
  agreement.}

  \item{loa.line.col}{color of the line for estimated
  limits of agreement.}

  \item{loa.text.col}{color of the label for estimated
  limits of agreement (defaults to the same as
  \code{loa.line.col}.)}

  \item{add.grid}{logical value. If \code{add.grid=TRUE}
  (Default) gridlines will be drawn.}

  \item{ylim}{limits for the y-axis}

  \item{cex}{numeric value specifying the magnification
  factor used for points}

  \item{...}{further graphical parameters}
}
\description{
Draw different Bland-Altman plot modifications (see
parameter \code{plot.type}).
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

    plotDifference( model ) # Default plot.type=3
    plotDifference( model, plot.type=5)
    plotDifference( model, plot.type=7, ref.line.lty=3, ref.line.col="green3" )
}
\references{
Bland, J. M., Altman, D. G. (1986) Statistical methods for
assessing agreement between two methods of clinical
measurement.  \emph{Lancet}, \bold{i:} 307--310.
}
\seealso{
\code{\link{plot.mcr}}, \code{\link{plotResiduals}},
\code{\link{plotDifference}}, \code{\link{plotBias}},
\code{\link{compareFit}}
}

