\name{MCResultJackknife-class}
\Rdversion{1.1}
\docType{class}
\alias{MCResultJackknife-class}
\alias{calcResponse,MCResultJackknife-method}
\alias{getRJIF,MCResultJackknife-method}
\alias{plotwithRJIF,MCResultJackknife-method}
\alias{printSummary,MCResultJackknife-method}

\title{Class \code{"MCResultJackknife"}}
\description{
Result of a method comparison with Jackknife based confidence intervals.
}
\section{Objects from the Class}{
Object is typically created by a call to function \code{\link{mcreg}}.
Object can be directly constructed by calling \code{\link{newMCResultJackknife}} or \code{new("MCResultJackknife", data, para, mnames, regmeth, cimeth, alpha, glob.coef, B0jack, B1jack, error.ratio, weight)}.
}
\section{Slots}{
  \describe{
    \item{\code{glob.coef}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{B0jack}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{B1jack}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{data}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{para}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{mnames}:}{Object of class \code{"character"} ~~ }
    \item{\code{regmeth}:}{Object of class \code{"character"} ~~ }
    \item{\code{cimeth}:}{Object of class \code{"character"} ~~ }
    \item{\code{error.ratio}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{alpha}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{weight}:}{Object of class \code{"numeric"} ~~ }
  }
}

\section{Extends}{
Class \code{"\linkS4class{MCResult}"}, directly.
}
\section{Methods}{
  \describe{
    \item{calcResponse}{\code{signature(.Object = "MCResultJackknife")}: ... }
    \item{getRJIF}{\code{signature(.Object = "MCResultJackknife")}: ... }
    \item{plotwithRJIF}{\code{signature(.Object = "MCResultJackknife")}: ... }
    \item{printSummary}{\code{signature(.Object = "MCResultJackknife")}: ... }
	 }
}
\author{
Ekaterina Manuilova \email{ekaterina.manuilova@roche.com}, 
Andre Schuetzenmeister \email{andre.schuetzenmeister@roche.com}, 
Fabian Model \email{fabian.model@roche.com}
}
\keyword{classes}
\examples{
showClass("MCResultJackknife")
}
