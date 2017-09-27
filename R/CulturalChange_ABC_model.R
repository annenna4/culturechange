
culturalChange_ABC_model<-function(theta, samples, timePoints, mutationRate, 
                     sMean, sVariance, alphaDirichlet, rMean = NULL, rVariance = NULL, popSize = NULL, n_step_for_sample)
#Output: generative model as a function of theta which serves as the input to the ABC procedure
{ 
 
  nSamples = length(samples) #Number of samples
  nType = length(samples[[1]])+1 #Number of variant types
  iniSample = c(samples[[1]],0) #vector of the initial sample + '0' for mutations
  
  sampleSize = rep(0,nSamples) #A vector containing the size of each sample
  for (i in 1:nSamples){
    sampleSize[i] = sum(samples[[i]])
  }
  
  if(length(timePoints) != nSamples){stop("There should be as many dates as samples!")}
  
  
  tSamples<-function(theta){
  
    #Generation of population size for each time step between t1 and t2
    if(is.null(popSize)){
      popSize = generate_popSize(sampleSize, nSamples, timePoints, sMean, sVariance)
    }
    if(length(popSize) != timePoints[nSamples]-timePoints[1]+1){
      print("error : popSize does not have the right size")
    }
     
    #Generation of population structure at t1 
    iniPop = round(generate_initialPop(iniSample, 1, alphaDirichlet) * popSize[1])
    popSize[1] = sum(iniPop)
    
    #Definition of b (strength of frequency-dependent selection) and r (fraction of population to be replaced in each time step)
    if(length(theta) == 2){ #Inference of both b and r from the data
      r=theta[2]
      b=theta[1]
    } 
    if(length(theta) == 1){ #Inference of b only
      r = rnorm(1,rMean,rVariance)
      while(r <= 0){
        r = rnorm(1,rMean,rVariance)
      }
      b = theta
    }

    #Generation of variants to be removed and added  
    h = generate_removalReplacement(popSize, timePoints, r, nSamples)
    u = h[[1]]
    v = h[[2]]  
  
    #Generation of theoretical samples at t1,t2,...,tn under transmission process specified by b
    theoSamples = generate_cultural_change(iniPop, u, v, b, mutationRate, timePoints, nType, sampleSize, nSamples, n_step_for_sample)
  
    return(theoSamples)
  }
  
  return(tSamples)
}

