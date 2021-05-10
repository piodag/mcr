\name{MCResult.calcCUSUM}
\alias{MCResult.calcCUSUM}
\alias{calcCUSUM}
\title{Calculate CUSUM Statistics According to Passing & Bablok (1983)}
\usage{
MCResult.calcCUSUM(.Object)
}
\arguments{
  \item{.Object}{object of class "MCResult".}
}
\value{
A list containing the following elements: \item{nPos}{sum
of positive residuals} \item{nNeg}{sum of negative
residuals} \item{cusum}{a cumulative sum of vector with
scores ri for each point, sorted increasing by distance of
points to regression line.} \item{max.cumsum}{Test
statisics of linearity test}
}
\description{
Calculate CUSUM Statistics According to Passing & Bablok
(1983)
}
\references{
Passing, H., Bablok, W. (1983) A new biometrical procedure
for testing the equality of measurements from two different
analytical methods.  Application of linear regression
procedures for method comparison studies in clinical
chemistry, Part I.  \emph{J Clin Chem Clin Biochem}.  Nov;
\bold{21(11)}:709--20.
}

