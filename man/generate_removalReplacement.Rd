
\name{generate_removalReplacement}
\alias{generate_removalReplacement}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculation of magnitude of population replacement in each time step}
\description{
 Generates the number of variants to be removed and to be added in each timestep.
} 
\usage{
generate_removalReplacement <- function(popSize, timePoints, r, nSamples)
}
\arguments{
    \item{popSize}{vector containing the the population sizes for each time step between t1 and tnSamples}
    \item{timePoints}{vector specifying the time points tj, j=1,...,nSamples}
    \item{r}{fraction of the population to be removed in each time step}
    \item{nSamples}{number of observed samples}
}

\details{
This function calculates the number of variants to be removed from the population (determined by fraction r) and the number of variants to be added so that the desired popualtion size is reached in each time step. 
}
%- maybe also 'usage' for other objects documented here.
\references{}
\author{Clementine Straub}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }