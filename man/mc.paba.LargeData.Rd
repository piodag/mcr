\name{mc.paba.LargeData}
\alias{mc.paba.LargeData}
\title{Passing-Bablok Regression for Large Datasets}
\usage{
mc.paba.LargeData(X, Y, NBins = 1e+06, alpha = 0.05, posCor = TRUE,
  calcCI = TRUE)
}
\arguments{
  \item{X}{(numeric) vector containing measurement values
  of reference method}

  \item{Y}{(numeric) vector containing measurement values
  of test method}

  \item{NBins}{(integer) value specifying the number of
  bins used to classify slope-values}

  \item{alpha}{(numeric) value specifying the
  100(1-alpha)\% confidence level for confidence intervals}

  \item{posCor}{(logical) should algorithm assume positive
  correlation, i.e. symmetry around slope 1?}

  \item{calcCI}{(logical) should confidence intervals be
  computed?}
}
\value{
Matrix of estimates and confidence intervals for intercept
and slope. No standard errors provided by this algorithm.
}
\description{
This function represents an interface to a fast
C-implementation of an adaption of the Passing-Bablok
algorithm for large datasets. Instead of building the
complete matrix of pair-wise slope values, a pre-defined
binning of slope-values is used (Default NBins=1e06). This
reduces the required memory dramatically and speeds up the
computation.
}
\examples{
library("mcr")
 data(creatinine,package="mcr")

# remove any NAs
crea <- na.omit(creatinine)

# call the approximative Passing-Bablok algorithm (Default NBins=1e06)
res1 <- mcreg(x=crea[,1], y=crea[,2], method.reg="PaBaLarge", method.ci="analytical")
getCoefficients(res1)

# now increase the number of bins and see whether this makes a difference
res2 <- mcreg(x=crea[,1], y=crea[,2], method.reg="PaBaLarge", method.ci="analytical", NBins=1e07)
getCoefficients(res2)
getCoefficients(res1)-getCoefficients(res2)
}
\author{
Andre Schuetzenmeister
\email{andre.schuetzenmeister@roche.com} (partly re-using
code of function 'mc.paba')
}

