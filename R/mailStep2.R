#' mailStep2
#' 
#' @export
#' 

mailStep2 <- function(numSelectionIter,numModels,xExp,yExp,firstSOILWeightType,firstSOILPsi,verbose=FALSE) {
  p <- dim(xExp)[2]
  NExp <- dim(xExp)[1]
  
  allSOILScores <- rep(0,times=p) ##
  
  ## added bootstrap on 7/29/2024
  if (numSelectionIter > 1) {
    soilScoreMat <- matrix(NA,nrow=numSelectionIter,ncol=p)
    soilRes <- NULL
    for (i in 1:numSelectionIter) {
      if (verbose == TRUE) {
        print(sprintf("\tStep 2: Iteration %d",i))
      }
      
      tempInds <- sample(1:NExp,size=NExp,replace=TRUE)
      tempXExp <- xExp[tempInds,]
      tempYExp <- yExp[tempInds]
      if (firstSOILWeightType != "ARM") {
        tempSOILRes <- SOIL(x=tempXExp,y=tempYExp,n_bound = numModels,
                            weight_type=firstSOILWeightType,
                            psi=firstSOILPsi,family="gaussian",method="union")
      }
      else {
        tempSOILRes <- SOIL(x=tempXExp,y=tempYExp,
                            weight_type = "ARM",
                            psi=firstSOILPsi,family="gaussian",method="union",
                            n_train = ceiling(NExp/2)+4)
      }
      soilScoreMat[i,] <- as.numeric(tempSOILRes$importance)
      allSOILScores <- allSOILScores + as.numeric(tempSOILRes$importance)
    }
    soilUncertaintyVec <- 1 + apply(soilScoreMat,2,sd,na.rm=TRUE)
    allSOILScores <- allSOILScores / numSelectionIter    
  }
  else {
    soilScoreMat <- NULL
    if (firstSOILWeightType != "ARM") {
      soilRes <- SOIL(x=xExp,y=yExp,n_bound = numModels,
                      weight_type=firstSOILWeightType,
                      psi=firstSOILPsi,family="gaussian",method="union")
    }
    else {
      soilRes <- SOIL(x=xExp,y=yExp,
                      weight_type = "ARM",
                      psi=firstSOILPsi,family="gaussian",method="union",
                      n_train = ceiling(NExp/2)+4)
    }
    soilUncertaintyVec <- rep(1,p)
    allSOILScores <- allSOILScores + as.numeric(soilRes$importance)    
  }
  
  resList <- list(allSOILScores = allSOILScores,
                  soilUncertaintyVec = soilUncertaintyVec,
                  soilRes = soilRes,
                  soilScoreMat = soilScoreMat)
  return(resList)
  
}