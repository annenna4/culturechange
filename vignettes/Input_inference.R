##############
##  Data    ##
##############

#absolute frequencies
S1 = c(1, 2, 5, 10, 20, 40, 60)
#S2 and S3 obtained with b = 0 and r=0.1
S2 = c(0, 3, 8, 16, 11, 43, 55, 2)
S3 = c(4,  0,  1, 19, 10, 44, 62, 7)

#S2 and S3 obtained with b = 0.1 and r=0.1
# S2 = c(0, 5, 1, 12, 13, 35, 72, 0)
# S3 = c(0, 2, 0, 4, 9, 27, 105, 0)

# S2 and S3 obtained with b = -0.1 and r=0.1
# S2 = c(5, 6, 8, 12, 30, 32, 40, 5)
# S3 = c(12, 8, 9, 13, 27, 31, 28, 19)

samples = list(S1,S2,S3)
#time points defining the temporal differences between samples
t1 = 1
t2 = 50
t3 = 100
timePoints = c(t1,t2,t3)


##############################
# Parameters for simulations #
##############################

mutationRate = 0.001
# definition of the population size
popSize = NULL
# s is the ratio between the sample size and the population size, normally distributed for each sample
# with mean sMean and variance sVariance
sMean = 0.1
sVariance = 0.05
# r is the fraction of the population being replaced at each time step, normally distributed
# with mean rMean and variance rVariance
rMean = 0.3
rVariance = 0
# parameters of the Dirichlet distribution used to create the initial population.
alphaDirichlet = rep(1,length(S1)+1)
# number of steps used to the sampling. Ex : if n_step_for_sample = 2, the simulated sample for t2 is created by randomly
# sampling half of the variants from t2 and half from (t2 - 1).
n_step_for_sample = 1

################################
# Parameters for ABC inference #
################################
#ABCMethod = "Rejection"
ABCMethod = "Delmoral"

# prior distribution for b and r
prior_b = c("unif",-0.2,0.2)
prior_r = c("unif",0.05,0.5)
prior_param = list(prior_b, prior_r)
#prior_param = list(prior_b)


# tolerance between theoretical samples (S2,S3) and simulated samples.
toleranceTargetABC = 0.01
# Rejection rate
tolRejection = 0.001
#positive integer specifying the number of simulations with an error value below the
#tolerance_target
nbSimulBelowTar = 1000
# specific parameters for Delmoral algorithm
alphaDelmoral = 0.9
MDelmoral = 5
