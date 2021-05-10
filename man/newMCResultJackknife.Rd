\name{newMCResultJackknife}
\alias{newMCResultJackknife}
\title{MCResultJackknife Object Constructor with Matrix in Wide Format as Input}
\usage{
newMCResultJackknife(wdata, para, sample.names = NULL, method.names = NULL,
  regmeth = "Unknown", glob.coef, cimeth = "unknown", B0jack, B1jack,
  error.ratio = error.ratio, alpha = 0.05, weight = rep(1, nrow(wdata)))
}
\arguments{
  \item{wdata}{measurement data in matrix format. First
  column reference method (x), second column test method
  (y).}

  \item{para}{regression parameters in matrix form. Rows:
  Intercept, Slope. Cols: EST, SE, LCI, UCI.}

  \item{sample.names}{names of individual data points, e.g.
  barcodes of measured samples.}

  \item{method.names}{names of reference and test method.}

  \item{regmeth}{name of statistical method used for
  regression.}

  \item{glob.coef}{global coefficients}

  \item{cimeth}{name of statistical method used for
  computing confidence intervals.}

  \item{B0jack}{jackknife intercepts}

  \item{B1jack}{jeckknife slopes}

  \item{error.ratio}{ratio between standard deviation of
  reference and test method.}

  \item{alpha}{numeric value specifying the
  100(1-\code{alpha})\% confidence level of confidence
  intervals (Default is 0.05).}

  \item{weight}{numeric vector specifying the weights used
  for each point}
}
\value{
MCResult object containing regression results.
}
\description{
MCResultJackknife Object Constructor with Matrix in Wide
Format as Input
}

