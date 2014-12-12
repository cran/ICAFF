\name{ICAFF-package}
\alias{ICAFF-package}
\alias{ICAFF}
\docType{package}
\title{
Imperialist Competitive Algorithm
}
\description{

Imperialist Competitive Algorithm starts by generating a set of candidate random solutions in the search space of the optimization problem.

The generated random points are called the initial countries.
Countries in this algorithm are the counterpart of Chromosomes in GAs and Particles in Particle Swarm Optimization (PSO) and it is an array of values of a candidate solution of optimization problem.

The cost function of the optimization problem determines the power of each country.

Based on their power, some of the best initial countries (the countries with the least cost function value), become imperialists and start taking control of other countries (called colonies) and form the initial empires.

Two main operators of this algorithm are Assimilation and Revolution.
Assimilation makes the colonies of each empire get closer to the imperialist state in the space of socio political characteristics (optimization search space).

Revolution brings about sudden random changes in the position of some of the countries in the search space.
 
During assimilation and revolution a colony might reach a better position and has the chance to take the control of the entire empire and replace the current imperialist state of the empire.

Imperialistic Competition is another part of this algorithm. 
All the empires try to win this game and take possession of colonies of other empires.
 
In each step of the algorithm, based on their power, all the empires have a chance to take control of one or more of the colonies of the weakest empire.

Algorithm continues with the mentioned steps (Assimilation, Revolution, Competition) until a stop condition is satisfied.
}
\details{
\tabular{ll}{
Package: \tab ICAFF\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-12-12\cr
License: \tab  GPL(>=2)\cr
}
}
\author{
Farimah Houshmand , Farzad Eskandari Ph.D.

Maintainer: Farimah Houshmand <hoshmandcomputer@gmail.com>
}
\references{

Atashpaz-Gargari, E. and Lucas, C.(2007). "Imperialist Competitive Algorithm: An algorithm for optimization inspired by imperialistic competition". IEEE Congress on Evolutionary Computation, Vol. 7,pp. 4661-4666.


Houshmand, F. and Eskandari, F.(2014)."Introduction The Imperialist Competitive Algorithm in R ". Proceeding of First Joint Seminar on Image processing and Intelligent Systems.Institute of Giasoldin Jamshid Kashani, Theran, Iran.


Jasour, A., Atashpaz-Gargari, E. and Lucas, C.(2008)."Vehicle Fuzzy Controller Design Using Imperialist Competitive Algorithm". In Second First Iranian Joint Congress on Fuzzy and Intelligent Systems, Tehran, Iran.


Jain, T. and Nigam, M.(2009)."Synergy of evolutionary algorithm and socio-political process for global optimization",Expert Systems with Applications. 


Rajabioun, R.Atashpaz-Gargari, E.and Lucas, C.(2008)."Colonial Competitive Algorithm as a Tool for Nash Equilibrium Point Achievement". Lecture Notes in Computer Science (LNCS).Springer 680-695.


Sepehri Rad, H. and Lucas, C.(2008)."Application of Imperialistic Competition Algorithm in Recommender Systems". 13th Int'l CSI Computer Conference (CSICC), Kish Island, Iran.

}
\keyword{
Imperialist, Competitive, optimization
}
\examples{
##  Example 1

##  cost function:f=x^2 +y^2+z^2+...

cost <- function(x){return(sum(x^2))}

ICA(cost,nvar=7, ncountries = 80, nimp = 10, maxiter = 100,lb = -10, ub= 10, 
           beta = 2, P_revolve = 0.3, zeta = 0.02)

##  Example 2

cost<-function(x){return(x[1]*sin(4*x[1])+ 1.1*x[2]*sin(2*x[2]))}

ICA(cost,2, ncountries = 80, nimp = 10, maxiter = 100,lb = -10, ub= 10, 
           beta = 2, P_revolve = 0.3, zeta = 0.02)

##  Example 3 

cost<-function(x){
return(0.5+(((sin(sqrt((x[1]^2)+(x[2]^2)))^2)-0.5)/(1+(0.1*((x[1]^2)
+(x[2]^2))))))}

ICA(cost,2, ncountries = 80, nimp = 10, maxiter = 100,lb = -10, ub= 10, 
           beta = 2, P_revolve = 0.3, zeta = 0.02)

##  Example 4 

cost<-function(x){
return(((((x[1]^2)+(x[2]^2))^0.25)*sin(30*((((x[1]+0.5)^2)+(x[2]^2))^0.1)))+
 (abs(x[1]))+ (abs(x[2])))}

ICA(cost,2)

##  Example 5

cost<-function(x){
return((-exp(-0.2*sqrt((x[1]^2)+(x[2]^2))))+(3*(cos(2*x[1])+sin(2*x[2]))))}

ICA(cost,2)

##  Example 6

cost<- function(x){
return((-x[1]*sin(sqrt(abs(x[1]-x[2]-9))))-((x[2]+9)*sin(sqrt(abs(x[2]+
(0.5*x[1])+9)))))}

ICA(cost,2)

}