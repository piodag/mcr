\name{MCResultResampling.initialize}
\alias{MCResultResampling.initialize}
\title{Initialize Method for 'MCResultAnalytical' Objects.}
\usage{
MCResultResampling.initialize(.Object, data = data.frame(X = NA, Y = NA),
  para = matrix(NA, ncol = 4, nrow = 2), xmean = 0, mnames = c("unknown",
  "unknown"), regmeth = "unknown", cimeth = "unknown",
  bootcimeth = "unknown", alpha = 0.05, glob.coef = c(0, 0),
  rng.seed = as.numeric(NA), rng.kind = "unknown", glob.sigma = c(0, 0),
  nsamples = 0, nnested = 0, B0 = 0, B1 = 0, MX = 0, sigmaB0 = 0,
  sigmaB1 = 0, error.ratio = 0, weight = 1, robust.cov = "MCD")
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

  \item{B0}{resampling intercepts}

  \item{B1}{resampling slopes}

  \item{MX}{Numeric vector with point estimations of
  (weighted-)average of reference method values for each
  bootstrap sample.}

  \item{sigmaB0}{SD for 'B0'}

  \item{sigmaB1}{SD for 'B1'}

  \item{weight}{1 for each data point}

  \item{robust.cov}{Robust covariance method for Box & Ellipses plot}
}
\description{
Method initializes newly created objects of class
'MCResultAnalytical'.
}

