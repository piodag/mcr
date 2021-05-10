\name{MCResult.getFitted}
\alias{MCResult.getFitted}
\alias{getFitted}
\title{Get Fitted Values.}
\usage{
MCResult.getFitted(.Object)
}
\arguments{
  \item{.Object}{object of class "MCResult".}
}
\value{
fitted values as data frame.
}
\description{
This funcion computes fitted values for a
'MCResult'-object. Depending on the regression method and
the error ratio, a projection onto the regression line is
performed accordingly. For each point (x_i; y_i) i=1,...,n
the projected point(x_hat_i; y_hat_i) is computed.
}
\seealso{
\code{\link{plotResiduals}} \code{\link{getResiduals}}
}

