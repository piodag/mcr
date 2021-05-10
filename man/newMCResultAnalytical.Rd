\name{newMCResultAnalytical}
\alias{newMCResultAnalytical}
\title{MCResultAnalytical object constructor with matrix in wide format as input.}
\usage{
newMCResultAnalytical(wdata, para, xmean, sample.names = NULL,
  method.names = NULL, regmeth = "Unknown", cimeth = "analytical",
  error.ratio = error.ratio, alpha = 0.05, weight = rep(1, nrow(wdata)))
}
\arguments{
  \item{wdata}{Measurement data in matrix format. First
  column reference method (x), second column comparator
  method (y).}

  \item{para}{Regression parameters in matrix form. Rows:
  Intercept, Slope. Cols: EST, SE, LCI, UCI.}

  \item{xmean}{Global (weighted) mean of x-values.}

  \item{sample.names}{Names of individual data points, e.g.
  barcodes of measured samples.}

  \item{method.names}{Names of reference and comparator
  method.}

  \item{regmeth}{Name of statistical method used for
  regression.}

  \item{cimeth}{Name of statistical method used for
  computing confidence intervals.}

  \item{error.ratio}{Ratio between standard deviation of
  reference and comparator method.}

  \item{alpha}{1 - significance level for confidence
  intervals.}

  \item{weight}{numeric vector specifying the weights used
  for each point}
}
\value{
MCResultAnalytical object containing regression results.
}
\description{
MCResultAnalytical object constructor with matrix in wide
format as input.
}

