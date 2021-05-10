\docType{package}
\name{mcr-package}
\alias{mcr}
\alias{mcr-package}
\title{Method Comparison Regression}
\description{
Method Comparison Regression. This package provides
regression methods to quantify the relation between two
measurement methods. In particular it addresses regression
problems with errors in both variables and without repeated
measurements. The package provides implementations of
Deming regression, weighted Deming regression, and
Passing-Bablok regression following the CLSI EP09-A3
recommendations for analytical method comparison and bias
estimation using patient samples.
}
\details{
The main function for performing regression analysis is
\code{\link{mcreg}}. Various functions for summarizing and
plotting regression results are provided (see examples in
\code{\link{mcreg}}).

For user site testing (installation verification) please
use the test case suite provided with the package. The test
case suite can be run by sourcing the 'runalltests.R'
script in the 'unitTests' folder. It requires the XML and
Runit packages.

\tabular{ll}{ Package: \tab mcr\cr Type: \tab Package\cr
Version: \tab 1.2.1\cr Date: \tab 2014-02-11\cr License: \tab
GPL 3\cr LazyLoad: \tab yes\cr }
}
\author{
Ekaterina Manuilova <ekaterina.manuilova@roche.com> Andre
Schuetzenmeister <andre.schuetzenmeister@roche.com> Fabian
Model <fabian.model@roche.com>
}
\references{
CLSI EP09-A3 \url{http://www.clsi.org}
}
\keyword{package}

