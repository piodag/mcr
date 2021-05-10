\name{MCResultBCa-class}
\Rdversion{1.1}
\docType{class}
\alias{MCResultBCa-class}
\alias{calcResponse,MCResultBCa-method}
\alias{printSummary,MCResultBCa-method}

\title{Class \code{"MCResultBCa"}}
\description{
Result of a method comparison with BCa-bootstrap based confidence intervals.
}
\section{Objects from the Class}{
Object is typically created by a call to function \code{\link{mcreg}}.
Object can be directly constructed by calling \code{\link{newMCResultBCa}} or \code{new("MCResultBCa", data, para, xmean, mnames, regmeth, cimeth, bootcimeth, alpha, glob.coef, glob.sigma, nsamples, nnested, B0jack, B1jack, B0, B1, MX, rng.seed, rng.kind, sigmaB0, sigmaB1, error.ratio, weight)}.
}
\section{Slots}{
  \describe{
    \item{\code{glob.sigma}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{xmean}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nsamples}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nnested}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{B0}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{B1}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{sigmaB0}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{sigmaB1}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{MX}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{bootcimeth}:}{Object of class \code{"character"} ~~ }
    \item{\code{rng.seed}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{rng.kind}:}{Object of class \code{"character"} ~~ }
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
    \item{\code{robust.cov}:}{Object of class \code{"logical"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{MCResultJackknife}"}, directly.
Class \code{"\linkS4class{MCResult}"}, by class "MCResultJackknife", distance 2.
}
\section{Methods}{
  \describe{
    \item{calcResponse}{\code{signature(.Object = "MCResultBCa")}: ... }
    \item{printSummary}{\code{signature(.Object = "MCResultBCa")}: ... }
	 }
}
\author{
Ekaterina Manuilova \email{ekaterina.manuilova@roche.com},
Andre Schuetzenmeister \email{andre.schuetzenmeister@roche.com},
Fabian Model \email{fabian.model@roche.com}
}
\examples{
showClass("MCResultBCa")
}
\keyword{classes}
