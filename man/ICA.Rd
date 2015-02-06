\name{ICA}
\alias{ICA}
\title{
Finding a minimum value for the optimization variables of a cost function.
}
\description{
ICA is a function for optimization by Imperialist Competitive Algorithm.
}
\usage{
ICA(cost, nvar, ncountries = 80, nimp = 10, maxiter = 100,
           lb = rep(-10, nvar), ub = rep(10, nvar), 
           beta = 2, P_revolve = 0.3, zeta = 0.02, \dots)
}
\arguments{
  \item{cost}{A cost function to be minimized.  The function accept the
    parameter values as a numerical vector as its principal argument.
    Additional arguments may be specified through the \code{\dots}
    argument below.}
  \item{nvar}{Number of optimization variables of cost function}
  \item{ncountries}{Number of initial countries}
  \item{nimp}{Number of Initial Imperialists}
  \item{maxiter}{Maximum number of iterations allowed.}
  \item{lb}{Lower limit of the optimization region; a numeric vector of
    length \code{nvar}.  Will be recycled if necessary.}
  \item{ub}{Upper limit of the optimization region; a numeric vector of
    length \code{nvar}.  Will be recycled if necessary.}
  \item{beta}{Assimilation coefficient.}
  \item{P_revolve}{Revolution is the process in which the
    socio-political characteristics of a country change suddenly.} 
  \item{zeta}{Total Cost of Empire = Cost of Imperialist + Zeta *
    mean(Cost of All Colonies)}
  \item{\dots}{Additional arguments, if needed, for the function \code{cost}.}
}
\details{
To use this code, you should only need to prepare your cost function.
}
\value{
  An object of class \code{"ICA"}, a list with components:
  \item{call}{The call used.}
  \item{postion}{The vector of components for the position of the
    minimum value found.}
  \item{value}{The minimum value at the optimal position.}
  \item{nimp}{The remaining number of imperialists at the conclusion of
    the procedure.}
  \item{trace}{A 1-column matrix of successive minimum values found at
    each iteration of the major loop.}
  \item{time}{The execution time taken to find the best solution.}
}
\references{
Atashpaz-Gargari, E. and Lucas, C. (2007). \emph{Imperialist Competitive
  Algorithm: An algorithm for optimization inspired by imperialistic
  competition}. IEEE Congress on Evolutionary Computation, Vol. 7,
  pp. 4661-4666. 
}
\author{Farimah Houshmand, Farzad Eskandari Ph.D. <Askandari@atu.ac.ir>

Maintainer: Farimah Houshmand <hoshmandcomputer@gmail.com>}

\note{
The steps can be summarized as the below pseudocode:

  0) Define objective function or cost function:  f(x, ...), x = (x[1], x[2],
  \dots, x[nvar]) ; 

  1) Initialization of the algorithm. Generate some random solution in
  the search space and create initial empires. 

  2) Assimilation: Colonies move towards imperialist states in different
  in directions. 

  3) Revolution: Random changes occur in the characteristics of some
  countries. 

  4) Position exchange between a colony and imperialist. A colony with a
  better position than the imperialist, has the chance to take the
  control of empire by replacing the existing imperialist. 

  5) Imperialistic competition: All imperialists compete to take
  possession of colonies of each other. 

  6) Eliminate the powerless empires. Weak empires lose their power
  gradually and they will finally be eliminated. 

  7) If the stop condition is satisfied, stop, if not go to 2.

  8) End

Assimilation:  Movement of colonies toward imperialists (Assimilation Policy)
Revolution:  A sudden change in the socio-political characteristics.
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
\keyword{optimize} % __ONLY ONE__ keyword per line
