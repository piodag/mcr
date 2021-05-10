\name{MCResultBCa.plotBoxEllipses}
\alias{MCResultBCa.plotBoxEllipses}
\title{Plot distriblution of bootstrap coefficients and Box and Ellipses for the hypotesis testing}
\usage{
MCResultBCa.plotBoxEllipses(.Object, robust.cov = c("MCD","SDe","Classic"))
}
\arguments{
  \item{.Object}{Object of class "MCResultBCa"}

  \item{robust.cov}{Robust covariance matrix method for ellipses.
  There are three options:\cr \code{"MCD"} - MCD robust covariance
  matrix.\cr \code{"SDe"} - MCD robust covariance matrix.\cr \code{"Classic"} - Classical
  non robust covariance matrix.\cr}
}
\description{
Plot distriblution of bootstrap coefficients (slope and
intercept), the covariance robust ellipses (Donoho - Stahel method)
and the box of the classical confidence intervals.
}

