library("parallel")

####################
# Predictive Check #
####################

# The two functions above perform the predictive check. 
# The first one performs cultural change from an initial population structure.
# The second one creates a matrix of several initial population structures and call the first function for each (with function parApply : simulations can be done with several cores)

cult_change_for_predictiveCheck<-function(population_structure)
# Input : a vector containing a population structure given by dirichlet distribution
# Output : a vector containing samples at t2...tn
{     
  # Generation of the initial population
  #Generation of population size for each time step between t1 and t2
  if(is.null(popSize)){
    popSize = generate_popSize(sampleSize, nSamples, timePoints, sMean, sVariance)
  }
  if(length(popSize) != timePoints[nSamples]-timePoints[1]+1){
    print("error : popSize does not have the right size")
  }
  initialPopulation = round(population_structure*popSize[1])
  
  #Definition of b and r
  b = sample(bCheck,1) #sample one b from the vector b
  if(is.null(rMean) == FALSE){ #if rMean is known, r is chosen according to the normal distribution
    r = rnorm(1,rMean,rVariance)
    while(r <= 0){
      r = rnorm(1,rMean,rVariance)
    }
  }else{ # if rMean is not known, r has been inferred with b
    r = sample(rCheck,1) #sample one r from the vector r
  }
  
  #Generation of the number of variants to be removed (u) and added (v)
  h = generate_removalReplacement(popSize, timePoints, r, nSamples)
  u = h[[1]]
  v = h[[2]]  
  
  #Generation of theoretical samples at t1,t2,...,tn under transmission process specified by b
  theoSamples = generate_cultural_change(initialPopulation, u, v, b, mutationRate, timePoints, nType, sampleSize, nSamples, n_step_for_sample)
  
  return(theoSamples)
}




predictiveCheck<-function(n_param, theta, samples,timePoints,mutationRate,popSize=NULL, n_step_for_sample, n_target_sim, n_cores)
# Output : a matrix containing n_target_sim simulations of cultural change
{
  #Define bTheo and rTheo
  if(dim(theta)[2] == 2){
    bCheck = theta[,2]
    rCheck = theta[,3]
  }else{
    bCheck = theta[,2]
    rCheck = NULL
  }

  
  nSamples = length(samples) #Number of samples
  nType = length(samples[[1]]) + 1 #Number of variant types
  initialSample = c(samples[[1]],0) #vector of the initial sample + '0' for mutations
  sampleSize = rep(0,nSamples) #A vector containing the size of each sample
  for (i in 1:nSamples){
    sampleSize[i] = sum(samples[[i]])
  }

  
  # Generation of the matrix containing 'n_target_sim' population structures
  population_structure=generate_initialPop(initialSample, n_target_sim, alphaDirichlet)
  
  
  # Needed for parallelisation
  if(exists("rMean") == FALSE){
    rMean = NULL    
    rVariance = NULL
  }
  
  
  # Parallelisation
  cl <- makeCluster(n_cores)  # Initiate cluster
  clusterExport(cl, list("sampleSize","bCheck","rCheck","nType","rMean","sMean","rVariance","sVariance",
                         "timePoints","nSamples","n_step_for_sample","mutationRate","popSize"), envir=environment())
  clusterExport(cl, list("generate_popSize","generate_cultural_change","generate_removalReplacement",
                         "generate_transmissionProb"))
  clusterEvalQ(cl, library(Matrix))
  clusterEvalQ(cl, library(gtools))
  
  samplesSIM = parApply(cl, population_structure, MARGIN = 1, FUN = cult_change_for_predictiveCheck) #perform simulation for each initial population and return the result in a matrix
  stopCluster(cl)
  
  
  # Norm and sort the matrix containing the simulated samples
  count = length(samples[[1]])+1
  samplesSIMN = samplesSIM
  #samplesSIMN=apply(samplesSIM,MARGIN=2,FUN=norme)
  for(i in 2:nSamples){
    samplesSIMN[seq((i - 2) * count + 1,(i - 2) * count + length(samples[[i]])),] = samplesSIM[seq((i - 2) * count + 1,(i - 2) * count + length(samples[[i]])),] / apply(samplesSIM[seq((i - 2) * count + 1,(i - 2) * count + length(samples[[i]])),],MARGIN = 2,FUN = sum)
  }
  samplesSIM_sorted = apply(samplesSIMN,MARGIN = 1,sort)
  
  #Calculate 95% prediction interval
  aux = round(n_target_sim * 0.025)
  if(aux==0){aux = 1}
  PIinf = samplesSIM_sorted[aux,]
  PIsup = samplesSIM_sorted[round(n_target_sim-n_target_sim*0.025),]
  IC=rbind(PIinf,PIsup) #
  
  maxi=max(PIsup) + 0.2

  # Display the 95% prediction interval
  par(mfrow = c(1,nSamples - 1))
  for(i in 2:nSamples){
    t = seq(1, length(samples[[i]]))
    plot(
      x = t,
      y = IC[2, seq((i - 2) * count + 1,(i - 2) * count + length(samples[[i]]))],
      pch = 15,
      col = 1,
      ylab = "Relative frequencies",
      xlab = "Variant types",
      ylim = c(0, maxi)
    )
    #lines(x = t, y = IC[2, seq((i - 2) * count + 1,(i - 2) * count + length(samples[[i]]))], col = 1)
    points(x = t, y = IC[1, seq((i - 2) * count + 1,(i - 2) * count + length(samples[[i]]))], pch = 15)
    #lines(x = t, y = IC[1, seq((i - 2) * count + 1,(i - 2) * count + length(samples[[i]]))], col = 1)
    
    sampleFreq = samples[[i]] / sum(samples[[i]])
    points(x = t,
           y = sampleFreq,
           pch = 20,
           col = 2)
    #lines(x = t, y = sampleFreq, col = 2)
    legend("topleft", legend = c("Prediction interval","Data"), pch = 16,col = seq(1,2))
  }
  title(main = "Prediction intervals (95%) \n of the marginal frequency distributions for each variant type",outer=TRUE,line=-2)
  
  l = list(matrix_simulations = samplesSIM,PredictionInterval = IC)
  return(l)
  
}


norme<-function(u)
{
  return(u / sum(u))
}


plot_posterior_predictive_check<- function(output, n_param, n_target_sim, n_cores,
                                           mutationRate, popSize, n_step_for_sample){
  intermediaryOutputs = read.table(output)
  if(n_param==2){
    #Display histograms for b and r
    par(mfrow = c(1,3))
    hb<-hist(intermediaryOutputs[,2],main="Posterior distribution of b",col = 4,xlab = "Strength of \n frequency-dependent selection")
    hr<-hist(intermediaryOutputs[,3],main="Posterior distribution of r",col = 3,xlab = "Fraction of the population \n to be replaced at each time step")
    plot(x=intermediaryOutputs[,2],y = intermediaryOutputs[,3],pch = 16,col = 1,
         xlab ='b', ylab ='r',main = "Joint posterior distribution")
    }else{
    #Display histogram for b
    par(mfrow = c(1,1))
    h<-hist(intermediaryOutputs[,2],main = "Posterior distribution of b",col = 4, xlab = "Strength of frequency-dependent selection")
    }
  predictiveCheckResults = predictiveCheck(n_param, intermediaryOutputs, samples,timePoints,
                                           mutationRate,popSize, n_step_for_sample, 
                                           n_target_sim, n_cores)
}
