\name{MCResultJackknife.initialize}
\alias{MCResultJackknife.initialize}
\title{Initialize Method for 'MCResultJackknife' Objects.}
\usage{
MCResultJackknife.initialize(.Object, data = data.frame(X = NA, Y = NA),
  para = matrix(NA, ncol = 4, nrow = 2), mnames = c("unknown", "unknown"),
  regmeth = "unknown", cimeth = "jackknife", alpha = 0.05,
  glob.coef = c(0, 0), B0jack = 0, B1jack = 0, error.ratio = 0,
  weight = 1)
}
\arguments{
  \item{.Object}{object to be initialized}

  \item{data}{empty data.frame}

  \item{para}{empty coefficient matrix}

  \item{mnames}{empty method names vector}

  \item{regmeth}{string specifying the regression-method}

  \item{cimeth}{string specifying the confidence interval
  method}

  \item{error.ratio}{for deming regression}

  \item{alpha}{value specifying the 100(1-alpha)%
  confidence-level}

  \item{glob.coef}{global coefficients}

  \item{B0jack}{jackknife intercepts}

  \item{B1jack}{jackknife slopes}

  \item{weight}{1 for each data point}
}
\description{
Method initializes newly created objects of class
'MCResultAnalytical'.
}

