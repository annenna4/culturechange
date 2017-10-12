
\name{generate_cultural_change}
\alias{generate_cultural_change}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of cultural change}
\description{
 Generates theoretical samples at times tj, j=2,...,nSamples under a specific process of cultural transmission.
} 
\usage{
generate_cultural_change<-function(iniPop, u, v, b, mutationRate, timePoints, nType, sampleSize, nSamples, n_step_for_sample)}
\arguments{
    \item{iniPop}{vector containing the relative population-level frequencies of the cultural variant types at t1}
    \item{u}{number of variants to be removed in every time step}
    \item{v}{number of variants to be added in every time step}
    \item{b}{strength of frequency-dependent selection}
    \item{mutationRate}{mutation rate per transmission event}
    \item{timePoints}{vector specifying the time points tj, j=1,...,nSamples}
    \item{nType}{number of variant types in the population at t1}
    \item{sampleSize}{vector of sample sizes at time points tj, j=1,...,nSamples}
    \item{nSamples}{number of observed samples}
    \item{n_step_for_sample}{number of (theoretical) populations (at time steps t(j-n_step_for_sample+1),...,tj) from which the samples at tj are evenly drawn}  
}

\details{
This functions generates the theoretical sample-level frequencies of all variant types present in t1 (and an additional type that contains variants not obserevd in the sample at t1) at times tj, j=2,...,nSamples under the cultural transmission process defined in function generate_transmissionProb.
}
%- maybe also 'usage' for other objects documented here.
\references{}
\author{Clementine Straub}
\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }