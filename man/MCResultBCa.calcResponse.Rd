\name{MCResultBCa.calcResponse}
\alias{MCResultBCa.calcResponse}
\title{Caluculate Response}
\usage{
MCResultBCa.calcResponse(.Object, x.levels, alpha = 0.05,
  bootcimeth = .Object@bootcimeth)
}
\arguments{
  \item{.Object}{object of class 'MCResultBCa'}

  \item{x.levels}{numeric vector specifying values of the
  reference method for which prediction should be made}

  \item{alpha}{significance level for confidence intervals}

  \item{bootcimeth}{character string specifying the method
  to be used for bootstrap confidence intervals}
}
\description{
Calculate predicted values for given values of the
reference-method.
}

