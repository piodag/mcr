\name{MCResult-class}
\Rdversion{1.1}
\docType{class}
\alias{MCResult-class}
\alias{calcBias,MCResult-method}
\alias{calcCUSUM,MCResult-method}
\alias{calcResponse,MCResult-method}
\alias{getCoefficients,MCResult-method}
\alias{getData,MCResult-method}
\alias{getErrorRatio,MCResult-method}
\alias{getRegmethod,MCResult-method}
\alias{getResiduals,MCResult-method}
\alias{getFitted,MCResult-method}
\alias{getWeights,MCResult-method}
\alias{plot,MCResult-method}
\alias{plotBias,MCResult-method}
\alias{plotDifference,MCResult-method}
\alias{plotResiduals,MCResult-method}
\alias{printSummary,MCResult-method}

\title{Class \code{"MCResult"}}
\description{
Result of a method comparison.
}
\section{Objects from the Class}{
Object is typically created by a call to function \code{\link{mcreg}}.
Object can be directly constructed by calling \code{\link{newMCResult}} or \code{new("MCResult", data, para, mnames, regmeth, cimeth, error.ratio, alpha, weight)}.
}
\section{Slots}{
  \describe{
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
\section{Methods}{
  \describe{
    \item{calcBias}{\code{signature(.Object = "MCResult")}: ... }
    \item{calcCUSUM}{\code{signature(.Object = "MCResult")}: ... }
    \item{calcResponse}{\code{signature(.Object = "MCResult")}: ... }
    \item{getCoefficients}{\code{signature(.Object = "MCResult")}: ... }
    \item{getData}{\code{signature(.Object = "MCResult")}: ... }
    \item{getErrorRatio}{\code{signature(.Object = "MCResult")}: ... }
    \item{getRegmethod}{\code{signature(.Object = "MCResult")}: ... }
    \item{getResiduals}{\code{signature(.Object = "MCResult")}: ... }
    \item{getWeights}{\code{signature(.Object = "MCResult")}: ... }
    \item{plot}{\code{signature(x = "MCResult")}: ... }
    \item{plotBias}{\code{signature(x = "MCResult")}: ... }
    \item{plotDifference}{\code{signature(.Object = "MCResult")}: ... }
    \item{plotResiduals}{\code{signature(.Object = "MCResult")}: ... }
    \item{printSummary}{\code{signature(.Object = "MCResult")}: ... }
	 }
}
\author{
Ekaterina Manuilova \email{ekaterina.manuilova@roche.com}, 
Andre Schuetzenmeister \email{andre.schuetzenmeister@roche.com}, 
Fabian Model \email{fabian.model@roche.com}
}
\examples{
showClass("MCResult")
}
\keyword{classes}

