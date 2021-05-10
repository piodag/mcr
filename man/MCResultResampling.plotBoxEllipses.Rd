\name{MCResultResampling.plotBoxEllipses}
\alias{MCResultResampling.plotBoxEllipses}
\title{Plot distriblution of bootstrap coefficients with robust covariance ellipses and a CI box}
\usage{
MCResultResampling.plotBoxEllipses(.Object, robust.cov = c("MCD","SDe","Classic"))
}
\arguments{
  \item{.Object}{Object of class "MCResultResampling"}

  \item{robust.cov}{Robust covariance matrix method for ellipses
  There are three options:\cr \code{"MCD"} - MCD robust covariance
  matrix.\cr \code{"SDe"} - MCD robust covariance matrix.\cr \code{"Classic"} - Classical
  non robust covariance matrix.\cr}
}
\description{
Plot distriblution of bootstrap coefficients (slope and
intercept) with robust covariance ellipses and a CI box fot
hypotesis testing.
}

