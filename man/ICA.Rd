\name{ICA}
\alias{ICA}
\title{
Finding a minimum value for the optimization variables of a cost function.
}
\description{
ICA is a function for optimization by Imperialist Competitive Algorithm.
}
\usage{
ICA(cost, nvar, ncountries = 80, nimp = 10, maxiter = 100,lb = -10, ub= 10, 
           beta = 2, P_revolve = 0.3, zeta = 0.02)
}
\arguments{
  \item{cost}{A function of cost}
  \item{nvar}{Number of optimization variables of your objective function}
  \item{ncountries}{Number of initial countries}
  \item{nimp}{Number of Initial Imperialists}
  \item{maxiter}{Number of iteration}
  \item{lb}{Lower limit of the optimization parameters}
  \item{ub}{Uper limit of the optimization parameters}
  \item{beta}{Assimilation coefficient is shown by "beta"}
  \item{P_revolve}{Revolution is the process in which the socio-political characteristics of a country change suddenly}
  \item{zeta}{Total Cost of Empire = Cost of Imperialist + Zeta * mean(Cost of All Colonies)}
}
\details{
To use this code, you should only prepare your cost function.
}
\value{
\item{Best Solution}{Best value for optimization variables}
\item{Best Fitness}{Minimum value for cost function}
\item{nimp}{Number of remained Imperialists}
}
\references{
Atashpaz-Gargari, E. and Lucas, C.(2007)."Imperialist Competitive Algorithm: An algorithm for optimization inspired by imperialistic competition". IEEE Congress on Evolutionary Computation, Vol. 7, pp. 4661-4666.
}
\author{Farimah Houshmand, Farzad Eskandari Ph.D.

Maintainer: Farimah Houshmand <hoshmandcomputer@gmail.com>}

\note{
The steps can be summarized as the below pseudocode:

  0) Define objective function or cost function:  f(x) , x= (x[1] , x[2] ,. . ., x[nvar]) ;

  1) Initialization of the algorithm. Generate some random solution in the search space and create initial empires.

  2) Assimilation: Colonies move towards imperialist states in different in directions.

  3) Revolution: Random changes occur in the characteristics of some countries.

  4) Position exchange between a colony and imperialist. A colony with a better position than the imperialist, has the chance to take the control of empire by replacing the existing imperialist.

  5) Imperialistic competition: All imperialists compete to take possession of colonies of each other.

  6) Eliminate the powerless empires. Weak empires lose their power gradually and they will finally be eliminated.

  7) If the stop condition is satisfied, stop, if not go to 2.

  8) End

Assimilation:  Movement of colonies toward imperialists (Assimilation Policy)
Revolution:  A sudden change in the socio-political characteristics.
}
\examples{

#--------cost function:f=x.sin(4x)+1.1y.sin(2y)

cost <- function(x){return(x[1]*sin(4*x[1])+ 1.1*x[2]*sin(2*x[2]))}

ICA(cost,nvar=2, ncountries = 80, nimp = 10, maxiter = 100,lb = -10, ub= 10, 
           beta = 2, P_revolve = 0.3, zeta = 0.02)
}
\keyword{ Assimilation,Revolution,ICA,Imperialist Competitive Algorithm,colony}
\keyword{ Imperialist Competitive Algorithm }% __ONLY ONE__ keyword per line
