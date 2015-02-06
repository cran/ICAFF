\name{plot.ICA}
\alias{plot.ICA}
\alias{print.ICA}
\alias{summary.ICA}
\alias{coef.ICA}
\title{
  Methods for ICA objects.
}
\description{
  Provide standard methods for manipulating ICA objects, namely
  printing, plotting, summarising and extracting the position of the minimum.
}
\usage{
\method{plot}{ICA}(x, \dots, xlab = "Iteration", ylab = "Value",
                   main = "ICA History", col = "red")
\method{print}{ICA}(x, \dots)
\method{summary}{ICA}(object, \dots)
\method{coef}{ICA}(object, \dots)
}
\arguments{
  \item{x, object}{An object of class \code{"ICA"}.}
  \item{xlab, ylab, main, col}{Graphics parameters.}
  \item{\dots}{Additional arguments passed on to the method.}
}
\details{
  Methods for standard generic functions when dealing with objects of
  class \code{"ICA"}
}
\value{
  \code{print} method: the value is printed and returned invisibly.

  \code{summary} method: dummy.  Returns the object unchanged.

  \code{plot} method: a plot of the history of the process is produced
  with a \code{NULL} return value.

  \code{coef} method: extract the location vector for the minimum value.
}

\examples{

## --------cost function: f(x,y) = x * sin(4 * x) + 1.1 * y * sin(2 * y)
## --------search region: -10 <= x, y <= 10

cost <- function(x) {
  x[1] * sin(4 * x[1]) + 1.1 * x[2] * sin(2 * x[2])
}

ICAout <- ICA(cost, nvar = 2, ncountries = 80, nimp = 10,
              maxiter = 100, lb = -10, ub = 10, 
              beta = 2, P_revolve = 0.3, zeta = 0.02)

summary(ICAout)     ## same as the print method
coef(ICAout)        ## get the position of the minimum
cost(coef(ICAout))  ## cost at the minimum
plot(ICAout)        ## show the history of the process

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize } % __ONLY ONE__ keyword per line
