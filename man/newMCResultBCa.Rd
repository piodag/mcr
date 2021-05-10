\name{newMCResultBCa}
\alias{newMCResultBCa}
\title{MCResultBCa object constructor with matrix in wide format as input.}
\usage{
newMCResultBCa(wdata, para, xmean, sample.names = NULL, method.names = NULL,
  regmeth = "unknown", glob.coef, glob.sigma, cimeth = "unknown",
  bootcimeth = "unknown", nsamples, nnested, rng.seed, rng.kind, B0jack,
  B1jack, B0, B1, MX, sigmaB0, sigmaB1, error.ratio, alpha = 0.05,
  weight = rep(1, nrow(wdata)))
}
\arguments{
  \item{wdata}{Measurement data in matrix format. First
  column reference method (x), second column comparator
  method (y).}

  \item{para}{Regression parameters in matrix form. Rows:
  Intercept, Slope. Cols: EST, SE, LCI, UCI.}

  \item{xmean}{Global (weighted) mean of x-values}

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

  \item{glob.coef}{Numeric vector of length two with global
  point estimations of intercept and slope.}

  \item{glob.sigma}{Numeric vector of length two with
  global estimations of standard errors of intercept and
  slope.}

  \item{bootcimeth}{Bootstrap based confidence interval
  estimation method.}

  \item{nsamples}{Number of bootstrap samples.}

  \item{nnested}{Number of nested bootstrap samples.}

  \item{rng.seed}{Seed used to call mcreg, NULL if no seed
  was used}

  \item{rng.kind}{RNG type (string, see set.seed for
  details) used, only meaningfull if rng.seed was
  specified}

  \item{B0jack}{Numeric vector with point estimations of
  intercept for jackknife samples.}

  \item{B1jack}{Numeric vector with point estimations of
  slope for jackknife samples.}

  \item{B0}{Numeric vector with point estimations of
  intercept for each bootstrap sample.}

  \item{B1}{Numeric vector with point estimations of slope
  for each bootstrap sample.}

  \item{sigmaB0}{Numeric vector with estimation of standard
  error of intercept for each bootstrap sample.}

  \item{sigmaB1}{Numeric vector with estimation of standard
  error of slope for each bootstrap sample.}

  \item{MX}{Numeric vector with point estimations of
  (weighted-)average of reference method values for each
  bootstrap sample.}

  \item{weight}{numeric vector specifying the weights used
  for each point}
}
\value{
MCResult object containing regression results.
}
\description{
MCResultBCa object constructor with matrix in wide format
as input.
}

