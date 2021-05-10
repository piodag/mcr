\name{MCResultBCa.initialize}
\alias{MCResultBCa.initialize}
\title{Initialize Method for 'MCResultBCa' Objects.}
\usage{
MCResultBCa.initialize(.Object, data = data.frame(X = NA, Y = NA),
  para = matrix(NA, ncol = 4, nrow = 2), xmean = 0, mnames = c("unknown",
  "unknown"), regmeth = "unknown", cimeth = "unknown",
  bootcimeth = "unknown", alpha = 0.05, glob.coef = c(0, 0),
  glob.sigma = c(0, 0), nsamples = 0, nnested = 0, B0jack = 0,
  B1jack = 0, B0 = 0, B1 = 0, MX = 0, rng.seed = as.numeric(NA),
  rng.kind = "unknown", sigmaB0 = 0, sigmaB1 = 0, error.ratio = 0,
  weight = 1, robust.cov = "MCD")
}
\arguments{
  \item{.Object}{object to be initialized}

  \item{data}{empty data.frame}

  \item{xmean}{0 for init-purpose}

  \item{para}{empty coefficient matrix}

  \item{mnames}{empty method names vector}

  \item{regmeth}{string specifying the regression-method}

  \item{cimeth}{string specifying the confidence interval
  method}

  \item{bootcimeth}{string specifying the method for
  bootstrap confidence intervals}

  \item{error.ratio}{for deming regression}

  \item{alpha}{value specifying the 100(1-alpha)%
  confidence-level}

  \item{glob.coef}{global coefficients}

  \item{rng.seed}{random number generator seed}

  \item{rng.kind}{type of the random number generator}

  \item{glob.sigma}{global sd values for regression
  parameters}

  \item{nsamples}{number of samples for resampling}

  \item{nnested}{number of inner simulation for nested
  bootstrap}

  \item{B0jack}{jackknife intercpet}

  \item{B1jack}{jackknife slope}

  \item{B0}{intercept}

  \item{B1}{slope}

  \item{MX}{parameter}

  \item{sigmaB0}{SD for intercepts}

  \item{sigmaB1}{SD for slopes}

  \item{weight}{1 for each data point}

  \item{robust.cov}{Robust covariance method for Box & Ellipses plot}
}
\description{
Method initializes newly created objects of class
'MCResultBCa'.
}

