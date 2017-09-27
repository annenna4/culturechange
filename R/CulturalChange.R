
generate_popSize<-function(sampleSize, nSamples, timePoints, sMean, sVariance)
#Output: vector N of size (tn-t1+1) containing an estimate of the popualtion size at each 
#        intermediate time step
#Estimation assumes that sample size is a fraction of population size at ti, i=1,...,n
#N(ti)=1/s*n(ti) with s being normally distributed with mean sMean and variance sVariance
#Population sizes at intermediate time points are obtained by linear interpolation
{
  #generating population sizes at ti 
  Nh = rep(0,nSamples) #vector containing the size of the population for each time point ti
  s = rep(0,nSamples) #vector containing the ratio between sample size and popluation size for each time point ti
  for (i in 1:nSamples){
    s[i] = rnorm(1,sMean,sVariance)
    while(s[i]<0){s[i]=rnorm(1,sMean,sVariance)}
    Nh[i] = round(sampleSize[i]/s[i])
  }
  #generating population sizes intermediate time points 
  N = rep(0, (timePoints[nSamples]-timePoints[1])+1)
  N[1] = Nh[1]
  for (i in 1:(nSamples-1)){
    time_diff = timePoints[i+1]-timePoints[i]+1
    #Linear interpolation between N[i] and N[i+1]
    h = (approx(c(1:2),c(Nh[i],Nh[i+1]),n = time_diff))
    N[(timePoints[i]-timePoints[1]+2):(timePoints[i+1]-timePoints[1]+1)] = round(h$y[2:time_diff])
  }
  return(N)
}

generate_initialPop<-function(iniSample, n_target_sim, alphaDirichlet)
#Output: Vector of size FIXME+1 containing relative frequencies of the variants present at
#        t1 in the population (estimation independent of total population size)
#Estimation is based on Dirichlet distribution
#Sample has to contains absolute frequencies   
{
  if(is.null(alphaDirichlet)){
    stop("No alpha values defined!")
  } else{
    if (isTRUE(0 %in% alphaDirichlet)){
      stop("All alpha values have to be larger than 0")
    }else{
      if(length(iniSample)!=length(alphaDirichlet)){
        stop("Length of alphaDirichlet needs to be (length of S1) + 1.")
      }
    }
    exponents = iniSample + alphaDirichlet
  }
  return(rdirichlet(n_target_sim,exponents))
}


generate_transmissionProb<-function(nType, pop, b, mutationRate)
#Output: transmission probability under the assumption of frequency-dependent selection where the strength of selection 
#        is controlled by the parameter b. b=0 indicates unbiased transmission. 
{
    h = rep(0,nType)
    transProb = c(0,(nType))
    #cb = nnzero(pop) #relative majority
    #cb = 0.3
    x = as.vector(pop / sum(pop))

    #h = x + (b * cb * x - b) #alternative form of frequency-dependent selection 
    h = x^(1 + b)
    
    #checks needed for alternative form
    #if frequency is 0 than transmssion probability needs to stay 0
    #index = which(x %in% 0)
    #h[index] = rep(0,length(index))
    #"negative" probabilities -> 0
    #index = which(h < 0)
    #h[index] = rep(0,length(index))
    
    #Normalisation
    h = h / sum(h)  
	  transProb[1:nType-1] = h[1:nType-1] * (1 - mutationRate)
	  transProb[nType] = h[nType] + (sum(h[1:nType-1])) * mutationRate
	  
	  return(transProb)	
}

generate_cultural_change<-function(iniPop, u, v, b, mutationRate, timePoints, nType, sampleSize, nSamples, n_step_for_sample)
#Output: Population-level frequencies of variants present in iniPop at ti under transmission process 
#        defined in generate_transmissionProb
#First step: random removal of u[t] variants
#Second step: Addition of v[t] variants according to specified  transmission process   
{
  theoSamples = rep(0,(nSamples - 1) * nType)
  pop = iniPop
  index = 1
  for (i in 1:(nSamples - 1)){
    n_step = timePoints[i+1] - timePoints[i] 
    if (i==(nSamples-1)){
      n_step = n_step + 1
    }
    xh = (i - 1) * nType + 1
    yh = i * nType
    q = sampleSize[i+1]%/%n_step_for_sample
    mod = sampleSize[i+1]%%n_step_for_sample
    for (t in 1:n_step){
      relFreq = pop / sum(pop)
      #Calculation of transmission probabilities
      transProb = generate_transmissionProb(nType,pop,b,mutationRate)
      if (u[(i-1) * (timePoints[i] - timePoints[1]) + t] > sum(pop)){
        stop("Number of variants to be removed is larger than population size")
      }
      if(u[(i-1) * (timePoints[i] - timePoints[1]) + t] < sum(pop)){
        h = sample(1:nType, u[(i-1) * (timePoints[i] - timePoints[1]) + t], replace = TRUE, prob = relFreq) 
        h = hist(h,seq(0.5,nType+0.5),plot = FALSE)
        poph = pop - h$counts
        index = which(poph < 0) #check whether there are negative values
        sumMinus = -sum(poph[index])
        #Repeating the process so that poph contains nonnegative entries
        while (sumMinus > 0){
          poph[index] = rep(0,length(index))
          pop[index] = rep(0,length(index))
          relFreq = pop / sum(pop)
      	  h = sample(1:nType, sumMinus, replace = TRUE, prob = relFreq) #choose the variants to remove
      	  h = hist(h,seq(0.5,nType + 0.5),plot = FALSE)		
      	  poph = poph - h$counts
      	  index = which(poph < 0)
      	  sumMinus = -sum(poph[index])
        }
      }else{
        poph = rep(0,length(pop))
      }
      #Addition of v[t] variants according to probability transProb
      h = sample(1:nType, v[(i-1) * (timePoints[i] - timePoints[1]) + t], replace = TRUE, prob = transProb)
      h = hist(h,seq(0.5,nType + 0.5),plot = FALSE)
      pop = poph + h$counts
      
      #Taking a sample from populations between t = (n_step - n_step_for_sample +1 ) and t = n_step
      if(t >= n_step - n_step_for_sample + 1 && t <= n_step - mod){
        n_to_take_at_each_step = q
        relFreq = pop / sum(pop)
        #random sampling of n variants from pop
        h = sample(1:nType,n_to_take_at_each_step,replace = TRUE,prob = relFreq)
        #determine relative frequencies of variants in sample 
        h = hist(h,seq(0.5,nType + 0.5),plot = FALSE)
        sam = h$counts
        theoSamples[xh:yh] = theoSamples[xh:yh] + sam
      }else if(t >= n_step - n_step_for_sample +1){
        n_to_take_at_each_step = q +1
        relFreq = pop / sum(pop)
        h = sample(1:nType,n_to_take_at_each_step,replace = TRUE,prob = relFreq)
        h = hist(h,seq(0.5,nType + 0.5),plot = FALSE)
        sam = h$counts
        theoSamples[xh:yh] = theoSamples[xh:yh] + sam
      }
    }
    theoSamples[xh:yh] = theoSamples[xh:yh] / sum(theoSamples[xh:yh])
    }
  return(theoSamples)
}

generate_removalReplacement <- function(popSize, timePoints, r, nSamples)
#Output: number of variants to be removed (u) and added (v) to the population in each time step between t1 and tn
{
  u = rep(0,(timePoints[nSamples] - timePoints[1]) + 1)
  v = rep(0,(timePoints[nSamples] - timePoints[1]) + 1)
  for (j in 1:(nSamples - 1)){
    n_step = timePoints[j+1] - timePoints[j] 
    if (j==(nSamples-1)){
      n_step = n_step 
    }
    for (i in 1:n_step) {
      if (popSize[(timePoints[j] - timePoints[1]) + i] <= popSize[(timePoints[j] - timePoints[1]) + i + 1]) {
        u[(timePoints[j] - timePoints[1]) + i] = floor(r * popSize[(timePoints[j] - timePoints[1]) + i])
        v[(timePoints[j] - timePoints[1]) + i] = popSize[(timePoints[j] - timePoints[1]) + i + 1] - (popSize[(timePoints[j] - timePoints[1]) + i] - u[(timePoints[j] - timePoints[1]) + i])
      } else{
        u[(timePoints[j] - timePoints[1]) + i] = floor(r * popSize[(timePoints[j] - timePoints[1]) + i + 1]) + (popSize[(timePoints[j] - timePoints[1]) + i] - popSize[(timePoints[j] - timePoints[1]) + i + 1])
        v[(timePoints[j] - timePoints[1]) + i] = popSize[(timePoints[j] - timePoints[1]) + i + 1] - (popSize[(timePoints[j] - timePoints[1]) + i] - u[(timePoints[j] - timePoints[1]) + i])
      }
    }
  }

  h = list(u,v)
  return(h)
}

