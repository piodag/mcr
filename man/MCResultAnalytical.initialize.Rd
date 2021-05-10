\name{MCResultAnalytical.initialize}
\alias{MCResultAnalytical.initialize}
\title{Initialize Method for 'MCResultAnalytical' Objects.}
\usage{
MCResultAnalytical.initialize(.Object, data = data.frame(X = NA, Y = NA),
  xmean = 0, para = matrix(NA, ncol = 4, nrow = 2), mnames = c("unknown",
  "unknown"), regmeth = "unknown", cimeth = "analytical", error.ratio = 0,
  alpha = 0.05, weight = 1)
}
\arguments{
  \item{.Object}{object to be initialized}

  \item{data}{empty data.frame}

  \item{xmean}{mean value}

  \item{para}{empty coefficient matrix}

  \item{mnames}{empty method names vector}

  \item{regmeth}{string specifying the regression-method}

  \item{cimeth}{string specifying the confidence interval
  method}

  \item{error.ratio}{for deming regression}

  \item{alpha}{value specifying the 100(1-alpha)%
  confidence-level}

  \item{weight}{1 for each data point}
}
\description{
Initialize Method for 'MCResultAnalytical' Objects.
}

