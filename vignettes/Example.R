remove (list=ls()) #clear workspace
source("CulturalChange.R")
source("culturalChange_ABC_model.R")
source("PredictiveCheck.R")
source("Input_inference.R")

############################
# Model of cultural change #
############################

# tdata is a function of theta. 
# tdata(theta) performs cultural change and returns the simulated samples at t2, t3...
tdata = culturalChange_ABC_model(theta, samples, timePoints, mutationRate, sMean, sVariance, 
                    alphaDirichlet, rMean, rVariance, popSize, n_step_for_sample)


#############################################
#### ABC analysis : infer theta #############
#############################################
empiricalSamples=NULL
for(i in 2:length(samples)){
empiricalSamples= c(empiricalSamples,unlist(samples[i])/sum(unlist(samples[i])))
}
n_param = length(prior_param)


if(ABCMethod == "Delmoral"){
    
    ABCResults = ABC_sequential(method="Delmoral", model = tdata, prior = prior_param,
                        nb_simul = nbSimulBelowTar, summary_stat_target = empiricalSamples,
                        tolerance_target = toleranceTargetABC, progress_bar=TRUE, 
                        alpha=alphaDelmoral, verbose=TRUE, M=MDelmoral)
    
    
    # ABC results
    plot_posterior_predictive_check("output_step2", n_param, n_target_sim=100, n_cores=2,
                                    mutationRate, popSize, n_step_for_sample)
    
    #Save the prediction interval in 'prediction_interval'
    write.table(format(t(predictiveCheckResults$PredictionInterval),digits=8),"prediction_interval.txt",
                sep="\t",row.name=FALSE,quote=FALSE)

}


if(ABCMethod == "Rejection"){
    ABCResultsRejection = ABC_rejection(model=tdata, prior = prior_param,
                                          nb_simul = nbSimulBelowTar,
                                          summary_stat_target=empiricalSamples,
                                          tol=tolRejection,verbose = TRUE, progress_bar = TRUE)
    
    nParamInferred = length(prior_param)
    #Cross Validation
    if(nParamInferred == 2){
      par(mfrow=c(2,1))
    }else{
      par(mfrow=c(1,1))
    }
    cv.rej <- cv4abc(param=ABCResultsRejection$param, sumstat=ABCResultsRejection$stats, nval=10,
                     tols=c(.1,.2,.3), method="rejection")
    plot(cv.rej)
    
    #Coverage
    resultRejection = read.table("output")
    coverageRes = cov.pi(as.data.frame(resultRejection[,1:nParamInferred]), as.data.frame(resultRejection[,-(1:nParamInferred)]), 
                         testsets=seq(1,100,1), tol=0.1, multicore = FALSE, cores, method = "rejection")
    hist(coverageRes$raw[,4],main = "Coverage plot for b", xlab = "b")
}

