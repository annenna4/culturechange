
\name{generate_transmissionProb}
\alias{generate_transmissionProb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculation of cultural transmission probabilities}
\description{
  Calculates the probability that a variant type is chosen to be copied.
} 
\usage{
generate_transmissionProb<-function(nType, pop, b, mutationRate)
}
\arguments{
    \item{nType}{number of variant types present in the population at t1}
    \item{pop}{relative frequencies of the variant types in the population}
    \item{b}{strength of frequency-dependent selection}
    \item{mutationRate}{mutation rate per transmission event}
}
\details{
This function determines the probability that a variant type is chosen to be copied under the assumption of frequency-dependent selection. The strength of selection is controlled by the parameter b whith b=0 indicating unbiased transmission. 
}
%- maybe also 'usage' for other objects documented here.
\references{}
\author{Clementine Straub}
\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }