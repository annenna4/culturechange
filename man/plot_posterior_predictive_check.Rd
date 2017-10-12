
\name{plot_posterior_predictive_check}
\alias{plot_posterior_predictive_check}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Display and analysis of inference output}
\description{
  Plots the marginal and joint posterior distributions of teh inferred parameter and carries out teh posterior predictive check.
}
\usage{
plot_posterior_predictive_check<-function(output, n_param, n_target_sim, n_cores, mutationRate, popSize, n_step_for_sample)
}
\arguments{
    \item{output}{joint posterior distribution of inferred parameters in .txt format}
    \item{nParam}{number of parameter to be inferred}
    \item{n_target_sim}{number of simulations}
    \item{n_cores}{number of cores used}
    \item{mutationRate}{mutation rate per transmission event}
    \item{popSize}{vector containing the size of the population for each time step}
    \item{n_step_for_sample}{number of (theoretical) populations (at time steps t(j-n_step_for_sample+1),...,tj) from which the samples at tj are evenly drawn}  
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