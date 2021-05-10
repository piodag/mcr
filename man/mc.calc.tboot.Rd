\name{mc.calc.tboot}
\alias{mc.calc.tboot}
\title{Bootstrap-t Method for Calculation of Resampling Confidence Intervals}
\usage{
mc.calc.tboot(Xboot, Sboot, xhat, shat, alpha)
}
\arguments{
  \item{Xboot}{vector of point estimates for bootstrap
  sample.  The i-th element contains the point estimate for
  the i-th bootstrap sample.}

  \item{Sboot}{vector of standard deviations for each
  bootstrap sample.  It schould be estimated with any
  analytical method or nonparametric with nested
  bootstrap.}

  \item{xhat}{point estimate for the complete data set
  (scalar).}

  \item{shat}{estimate of standard deviation for the
  complete data set (scalar).}

  \item{alpha}{numeric value specifying the 100(1-alpha)\%
  confidence level for the confidence interval (Default is
  0.05).}
}
\value{
a list with elements \item{est}{point estimate for the
complete data set (xhat).} \item{se}{estimate of standard
deviation for the complete data set (shat).}
\item{CI}{confidence interval for the point estimate.}
}
\description{
Calculate resampling confidence intervals for intercept,
slope or bias with t-Boot method given a vector of
bootstrap point estimates and a vector of bootstrap
standard deviations.
}
\references{
Carpenter, J., Bithell, J. (2000) Bootstrap confidence
intervals: when, which, what? A practical guide for medical
statisticians.  \emph{Stat Med}, \bold{19 (9)}, 1141--1164.
}

