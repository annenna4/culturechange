
\name{predictiveCheck}
\alias{predictiveCheck}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of cultural change}
\description{
Generates theoretical samples at times tj, j=2,...,nSamples under specific processes of cultural transmission soecified by the joint posterior distribution.  
} 
\usage{
predictiveCheck<-function(n_param, theta, samples, timePoints, mutationRate, popSize = NULL, n_step_for_sample, n_target_sim, n_cores)
}
\arguments{
    \item{n_param}{number of paramters to be inferred}
    \item{theta}{joint posterior distribution of inferred parameters}
    \item{samples}{observed samples at times tj, j=1,...,tnSamples}
    \item{timePoints}{vector specifying the time points tj, j=1,...,nSamples}
    \item{mutationRate}{mutation rate per transmission event}
    \item{popSize}{vector containing the size of the population for each time step}
    \item{n_step_for_sample}{number of (theoretical) populations (at time steps t(j-n_step_for_sample+1),...,tj) from which the samples at tj are evenly drawn}  
    \item{n_target_sim}{number of simulations}
    \item{n_cores}{number of cores used}
}
\details{

}
%- maybe also 'usage' for other objects documented here.
\references{}
\author{Clementine Straub}
\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }