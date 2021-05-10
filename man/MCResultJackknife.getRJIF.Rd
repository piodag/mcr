\name{MCResultJackknife.getRJIF}
\alias{MCResultJackknife.getRJIF}
\alias{getRJIF}
\title{Relative Jackknife Influence Function}
\usage{
MCResultJackknife.getRJIF(.Object)
}
\arguments{
  \item{.Object}{object of class "MCResultJackknife" or
  "MCResultResampling".}
}
\value{
a list of the following elements: \item{slope}{numeric
vector containing the values of relative jackknife function
of slope.} \item{intercept}{numeric vector containing the
values of relative jackknife function of intercept.}
}
\description{
Calculate the value of relative jackknife function for each
observation.
}
\references{
Efron, B. (1990) Jackknife-After-Bootstrap Standard Errors
and Influence Functions.  Technical Report , \bold{N 134}.
}

