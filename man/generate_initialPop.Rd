
\name{generate_initialPop}
\alias{generate_initialPop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculation of population structure at t1}
\description{
  Calculates possible population structures at time t1 so that the observed sample could have been drawn at random from these structures. 
} 
\usage{
generate_initialPop<-function(iniSample, n_target_sim, alphaDirichlet)}
\arguments{
    \item{iniSample}{vector containing the absolute frequencies of the observed cultural variant types at t1}
    \item{n_target_sim}{number of population structures to be generated}
    \item{alphaDirichlet}{vector of size (k+1), where k stands for the number of variant types present at t1, containing the parameters of the Dirichlet distribution}
}
\details{
This function will generate a vector of size (k+1) containing relative frequencies of the variant types present at t1 in the population. The variant type (k+1) contains variants not observed in the sample. The estimation is based on the Dirichlet distribution and is independent of total population size. The resulting population structure has the property that the observed sample at t1 could have been drawn at random from it. Importantly, iniSample has to contain absolute frequencies.  
}
%- maybe also 'usage' for other objects documented here.
\references{}
\author{Clementine Straub}
\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }