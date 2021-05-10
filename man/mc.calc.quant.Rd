\name{mc.calc.quant}
\alias{mc.calc.quant}
\title{Quantile Calculation for BCa}
\usage{
mc.calc.quant(X, alpha)
}
\arguments{
  \item{X}{numeric vector.}

  \item{alpha}{probabilty}
}
\value{
alpha-quantile of vector X.
}
\description{
We are using the R default (SAS (type=3) seems bugged)
quantile calculation instead of the quantile function
described in Effron&Tibshirani.
}

