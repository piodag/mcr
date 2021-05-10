\name{MCResultAnalytical-class}
\Rdversion{1.1}
\docType{class}
\alias{MCResultAnalytical-class}
\alias{calcResponse,MCResultAnalytical-method}
\alias{printSummary,MCResultAnalytical-method}

\title{Class \code{"MCResultAnalytical"}}
\description{
Result of a method comparison based on analytical methods for computing confidence intervals.
}
\section{Objects from the Class}{
Object is typically created by a call to function \code{\link{mcreg}}.
Object can be directly constructed by calling \code{\link{newMCResultAnalytical}} or \code{new("MCResultAnalytical", data, xmean, para, mnames, regmeth, cimeth, error.ratio, alpha, weight)}.
}
\section{Slots}{
  \describe{
    \item{\code{xmean}:}{Object of class \code{"numeric"} ~~ }
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
    \item{calcResponse}{\code{signature(.Object = "MCResultAnalytical")}: ... }
    \item{printSummary}{\code{signature(.Object = "MCResultAnalytical")}: ... }
	 }
}
\author{
Ekaterina Manuilova \email{ekaterina.manuilova@roche.com}, 
Andre Schuetzenmeister \email{andre.schuetzenmeister@roche.com}, 
Fabian Model \email{fabian.model@roche.com}
}
\examples{
showClass("MCResultAnalytical")
}
\keyword{classes}
