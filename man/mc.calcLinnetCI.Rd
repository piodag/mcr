\name{mc.calcLinnetCI}
\alias{mc.calcLinnetCI}
\title{Jackknife Confidence Interval}
\usage{
mc.calcLinnetCI(Xjack, xhat, alpha = 0.05)
}
\arguments{
  \item{Xjack}{vector of point estimates for jackknife
  samples.  The i-th element contains point estimate for
  the dataset without the i-th observation.}

  \item{xhat}{point estimate for the complete data set
  (scalar).}

  \item{alpha}{numeric value specifying the 100(1-alpha)\%
  confidence level for the confidence interval (Default is
  0.05).}
}
\value{
a list with elements \item{est}{point estimate for the
complete data set (scalar).} \item{se}{standard deviation
of point estimate calculated with Jackknife Method.}
\item{CI}{confidence interval for point estimate.}
}
\description{
Calculate Jackknife confidence intervals for intercept,
slope or bias given of vector of jackknife point estimates
and global point estimate.
}
\references{
Linnet, K. (1993) Evaluation of Regression Procedures for
Methods Comparison Studies.  \emph{CLIN. CHEM.}
\bold{39/3}, 424--432.
}

