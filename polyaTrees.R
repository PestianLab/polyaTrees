#######################################
#
# polyaTrees 
#
# Credit:  Brian Connolly / Pestian Lab
# License:  MIT (https://github.com/PestianLab/polyaTrees/blob/master/LICENSE)
# Copyright (c) 2018-present, Pestian Lab / Cincinnati Children's Hospital Medical Center
#
#######################################

# based on Holmes, Chris C., et al. "Two-sample Bayesian nonparametric hypothesis testing." Bayesian Analysis 10.2 (2015): 297-320.
softmax<-function(x) { return(exp(x)/(1+exp(x))) }

# returns probability test distribution is a case
polyaTwoSampleComparison<-function(caseDistribution,controlDistribution,testDistribution,mapDistributionsToSoftmax,maxLevels=100000,prior=1,epsilon=0,verbose=FALSE) {
 if (epsilon>0) {
    caseDistribution<-caseDistribution+rnorm(length(caseDistribution),0,epsilon)
    controlDistribution<-controlDistribution+rnorm(length(controlDistribution),0,epsilon)
    testDistribution<-testDistribution+rnorm(1,0,epsilon)
 }
 # OK, so if mapping distribution to softmax, then distribution values between 0 and 1.  However, if I'm not mapping it, 
 # then I calculate the bounds of the distributions by increasing the boundaries by factors of two until it covers all the points of all the distributions -- hokey, but whatever.
 # start by 
 if (mapDistributionsToSoftmax) {
    totalDistribution<-c(caseDistribution,controlDistribution)
    if (testDistribution>max(totalDistribution)) { 
        return(1)
    } else if (testDistribution<min(totalDistribution)) {
        return(0)
    }
    minTotal<-min(totalDistribution)
    maxTotal<-max(totalDistribution)
    caseDistribution<-(caseDistribution-minTotal)/(maxTotal-minTotal)
    controlDistribution<-(controlDistribution-minTotal)/(maxTotal-minTotal)
    testDistribution<-(testDistribution-minTotal)/(maxTotal-minTotal)
    minX<-0;
    maxX<-1;
 } else {
    totalDistribution<-c(caseDistribution,controlDistribution)
    medianOfTotal<-median(c(caseDistribution,controlDistribution))
    minX<-medianOfTotal-max(max(totalDistribution)-medianOfTotal,medianOfTotal-min(totalDistribution))
    maxX<-medianOfTotal+max(max(totalDistribution)-medianOfTotal,medianOfTotal-min(totalDistribution))+0.000000000001
    if (verbose) print(paste("minX=",minX,"maxX=",maxX,"testDistribution=",testDistribution))
    if (verbose) print(paste(testDistribution,minX,maxX))
    if (testDistribution>maxX) {
       return(1)
    } else if (testDistribution<minX) {
       return(0)
    }

    # create polya tree levels - create list
    breaks<-c(minX,minX+(1:(2^1))*(maxX-minX)/(2^1))

    polyaLevelsCase <- list()
    polyaLevelsCase[[1]] <- hist(caseDistribution,breaks=breaks,plot=FALSE)$counts
    if (verbose) print("polyaLevelsCase[[1]]")
    if (verbose) print(polyaLevelsCase[[1]])

    polyaLevelsControl <- list()
    polyaLevelsControl[[1]] <- hist(controlDistribution,breaks=breaks,plot=FALSE)$counts

    polyaLevelsTest <- list()
    polyaLevelsTest[[1]] <- hist(testDistribution,breaks=breaks,plot=FALSE)$counts

    jLevel<-0;
    # alright.  So keep breaking up distributions until they are split into the number of unique values. 
    totalDistribution<-c(caseDistribution,controlDistribution,testDistribution)
    continue<-TRUE
    while ( continue ) {
        jLevel=jLevel+1
        breaks<-c(minX,minX+(1:(2^jLevel))*(maxX-minX)/(2^jLevel))
        polyaLevelsCase[[jLevel]] <- hist(caseDistribution,breaks=breaks,plot=FALSE)$counts
        polyaLevelsControl[[jLevel]] <- hist(controlDistribution,breaks=breaks,plot=FALSE)$counts
        polyaLevelsTest[[jLevel]] <- hist(testDistribution,breaks=breaks,plot=FALSE)$counts
        polyaLevelsTotal <- hist(totalDistribution,breaks=breaks,plot=FALSE)$counts
        if (verbose) print(paste(sum(polyaLevelsTotal>0),length(unique(totalDistribution))))
        continue<-(sum(polyaLevelsTotal>0)<length(unique(totalDistribution))) && jLevel<=maxLevels
    }
    nPolyaLevels<-jLevel-1

    # loop through all the pairs for all the levels
    logbCase<-0; # for hypothesis that A is similar to C
    logbControl<-0; # for hypothesis that B is similar to C
    logb<-0
    for (iLevel in 1:nPolyaLevels) {
       if (verbose) print(paste("evaluating level",iLevel,"of",nPolyaLevels))
       # loop over pairs
       if (prior==1) alpha<-iLevel^2 else alpha<-2^-iLevel
       for ( iPair in seq(1,length(polyaLevelsCase[[iLevel]]),2) ) {  
            logbCase<-logbCase+lgamma(alpha+polyaLevelsCase[[iLevel]][iPair]+polyaLevelsTest[[iLevel]][iPair])+lgamma(alpha+polyaLevelsCase[[iLevel]][iPair+1]+polyaLevelsTest[[iLevel]][iPair+1])
            logbCase<-logbCase-lgamma(alpha+polyaLevelsCase[[iLevel]][iPair]+polyaLevelsTest[[iLevel]][iPair]+alpha+polyaLevelsCase[[iLevel]][iPair+1]+polyaLevelsTest[[iLevel]][iPair+1])
            logbCase<-logbCase+lgamma(alpha+polyaLevelsControl[[iLevel]][iPair])+lgamma(alpha+polyaLevelsControl[[iLevel]][iPair+1])-lgamma(alpha+polyaLevelsControl[[iLevel]][iPair]+alpha+polyaLevelsControl[[iLevel]][iPair+1]) # B is hangin' out there

            logbControl<-logbControl+lgamma(alpha+polyaLevelsControl[[iLevel]][iPair]+polyaLevelsTest[[iLevel]][iPair])+lgamma(alpha+polyaLevelsControl[[iLevel]][iPair+1]+polyaLevelsTest[[iLevel]][iPair+1])
            logbControl<-logbControl-lgamma(alpha+polyaLevelsControl[[iLevel]][iPair]+polyaLevelsTest[[iLevel]][iPair]+alpha+polyaLevelsControl[[iLevel]][iPair+1]+polyaLevelsTest[[iLevel]][iPair+1])
            logbControl<-logbControl+lgamma(alpha+polyaLevelsCase[[iLevel]][iPair])+lgamma(alpha+polyaLevelsCase[[iLevel]][iPair+1])-lgamma(alpha+polyaLevelsCase[[iLevel]][iPair]+alpha+polyaLevelsCase[[iLevel]][iPair+1]) # A is hangin' out there
       }
    }
    if (verbose) print(paste(testDistribution,1/(1+exp(logbControl-logbCase+log(length(controlDistribution))-log(length(caseDistribution)) ))))
    return(1/(1+exp(logbControl-logbCase+log(length(controlDistribution))-log(length(caseDistribution)))))
  }
}
