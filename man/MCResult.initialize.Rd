\name{MCResult.initialize}
\alias{MCResult.initialize}
\title{MCResult Object Initialization}
\usage{
MCResult.initialize(.Object, data = data.frame(X = NA, Y = NA),
  para = matrix(NA, ncol = 4, nrow = 2), mnames = c("unknown", "unknown"),
  regmeth = "unknown", cimeth = "unknown", error.ratio = 0,
  alpha = 0.05, weight = 1)
}
\arguments{
  \item{.Object}{object of class "MCResult"}

  \item{data}{measurement data in matrix format. First
  column reference method (x), second column test method
  (y).}

  \item{para}{regression parameters in matrix form. Rows:
  Intercept, Slope. Cols: EST, SE, LCI, UCI.}

  \item{mnames}{names of reference and test method.}

  \item{regmeth}{name of statistical method used for
  regression.}

  \item{cimeth}{name of statistical method used for
  computing confidence intervals.}

  \item{error.ratio}{ratio between standard deviation of
  reference and test method.}

  \item{alpha}{numeric value specifying the
  100(1-\code{alpha})\% confidence level of confidence
  intervals (Default is 0.05).}

  \item{weight}{weights to be used for observations}
}
\value{
MCResult object with initialized parameter.
}
\description{
MCResult Object Initialization
}

