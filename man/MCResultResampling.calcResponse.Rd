\name{MCResultResampling.calcResponse}
\alias{MCResultResampling.calcResponse}
\title{Caluculate Response}
\usage{
MCResultResampling.calcResponse(.Object, x.levels, alpha = 0.05,
  bootcimeth = .Object@bootcimeth)
}
\arguments{
  \item{.Object}{object of class 'MCResultResampling'}

  \item{x.levels}{numeric vector specifying values of the
  reference method for which prediction should be made}

  \item{alpha}{significance level for confidence intervals}

  \item{bootcimeth}{bootstrap confidence interval method to
  be used}
}
\description{
Calculate predicted values for given values of the
reference-method.
}

