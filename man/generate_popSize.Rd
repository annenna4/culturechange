
\name{generate_popSize}
\alias{generate_popSize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calcualtion of population size}
\description{
  Calculates the size of the population for each time step based on the observed sample size. 
}
\usage{
generate_popSize<-function(sampleSize, nSamples, timePoints, sMean, sVariance)
}
\arguments{
    \item{sampleSize}{vector of sample sizes at time points tj, j=1,...,nSamples}
    \item{timePoints}{vector specifying the time points tj, j=1,...,nSamples}
    \item{nSamples}{number of observed samples}
    \item{sMean}{mean value of the ratio s between sample size and population size}
    \item{sVariance}{variance of the ratio s between sample size and population size}
}
\details{
    This function will generate a vector of size (tnSamples-t1+1) containing an estimate of the population size at each intermediate time step. The estimation assumes that the sample size is a fraction of population size at tj, j=1,...,nSamples. In detail it holds N(tj)=1/s*n(tj) with s being normally distributed with mean sMean and variance sVariance. Population sizes at intermediate time points are obtained by linear interpolation.
}
%- maybe also 'usage' for other objects documented here.
\references{}
\author{Clementine Straub}
\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }