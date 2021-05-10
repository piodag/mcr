\name{mc.bootstrap}
\alias{mc.bootstrap}
\title{Resampling estimation of regression parameters and standard  errors.}
\usage{
mc.bootstrap(method.reg = c("LinReg", "WLinReg", "Deming", "WDeming","MDeming","MMDeming", "PaBa",
  "PaBaLarge"), X, Y, error.ratio, nsamples = 1000, jackknife = TRUE,
  bootstrap = c("none", "bootstrap", "nestedbootstrap"), nnested = 25,
  iter.max = 30, threshold = 1e-08, NBins = 1e+06,
  slope.measure = c("radian", "tangent"))
}
\arguments{
  \item{X}{Measurement values of reference method}

  \item{Y}{Measurement values of test method}

  \item{error.ratio}{Ratio between squared measurement
  errors of reference- and test method, necessary for
  Deming regression. Default 1.}

  \item{method.reg}{Regression method.  It is possible to
  choose between five regression types: \code{"LinReg"} -
  ordinary least square regression, \code{"WLinReg"} -
  weighted ordinary least square regression,\code{"Deming"}
  - Deming regression, \code{"WDeming"} - weighted Deming
  regression, \code{"Meming"} - weighted M-Deming
  regression, \code{"MMDeming"} - weighted MM-Deming
  regression, \code{"PaBa"} - Passing-Bablok regression.}

  \item{bootstrap}{Bootstrap based confidence interval
  estimation method.}

  \item{jackknife}{Logical value. If TRUE - Jackknife based
  confidence interval estimation method.}

  \item{nsamples}{Number of bootstrap samples.}

  \item{nnested}{Number of nested bootstrap samples.}

  \item{iter.max}{maximum number of iterations for weighted
  Deming iterative algorithm.}

  \item{threshold}{Numerical tolerance for weighted Deming
  iterative algorithm convergence.}

  \item{NBins}{number of bins used when
  'reg.method="PaBaLarge"' to classify each slope in one of
  'NBins' bins of constant slope angle covering the range
  of all slopes.}

  \item{slope.measure}{angular measure of pairwise slopes
  used for exact PaBa regression (see \code{\link{mcreg}}
  for details).\cr  \code{"radian"} - for data sets with even
  sample numbers median slope is calculated as average of
  two central slope angles.\cr \code{"tangent"} - for data sets
  with even sample numbers median slope is calculated as
  average of two central slopes (tan(angle)).\cr}
}
\value{
a list consisting of \item{glob.coef}{Numeric vector of
length two with global point estimations of intercept and
slope.} \item{glob.sigma}{Numeric vector of length two with
global estimations of standard errors of intercept and
slope.} \item{xmean}{Global (weighted-)average of reference
method values.} \item{B0jack}{Numeric vector with point
estimations of intercept for jackknife samples.  The i-th
element contains point estimation for data set without i-th
observation} \item{B1jack}{Numeric vector with point
estimations of slope for jackknife samples.  The i-th
element contains point estimation for data set without i-th
observation} \item{B0}{Numeric vector with point
estimations of intercept for each bootstrap sample.  The
i-th element contains point estimation for i-th bootstrap
sample.} \item{B1}{Numeric vector with point estimations of
slope for each bootstrap sample.  The i-th element contains
point estimation for i-th bootstrap sample.}
\item{MX}{Numeric vector with point estimations of
(weighted-)average of reference method values for each
bootstrap sample.  The i-th element contains point
estimation for i-th bootstrap sample.}
\item{sigmaB0}{Numeric vector with estimation of standard
error of intercept for each bootstrap sample.  The i-th
element contains point estimation for i-th bootstrap
sample.} \item{sigmaB1}{Numeric vector with estimation of
standard error of slope for each bootstrap sample.  The
i-th element contains point estimation for i-th bootstrap
sample.} \item{nsamples}{Number of bootstrap samples.}
\item{nnested}{Number of nested bootstrap samples.}
\item{cimeth}{Method of confidence interval calculation
(bootstrap).} \item{npoints}{Number of observations.}
}
\description{
Generate jackknife or (nested-) bootstrap replicates of a
statistic applied to data. Only a nonparametric ballanced
design is possible. For each sample calculate point
estimations and standard errors for regression
coefficients.
}
\author{
Ekaterina Manuilova
\email{ekaterina.manuilova@roche.com}, Fabian
Model \email{fabian.model@roche.com}
}
\references{
Efron, B., Tibshirani, R.J. (1993) \emph{An Introduction to
the Bootstrap}. Chapman and Hall.  Carpenter, J., Bithell,
J. (2000) Bootstrap confidence intervals: when, which,
what? A practical guide for medical statisticians.
\emph{Stat Med}, \bold{19 (9)}, 1141--1164.
}

