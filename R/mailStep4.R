mailStep4 <- function(allSOILScores,selectedSet,xExp,yExp,smallestModelWeightType,smallestModelPsi,numModels) {
  selectedSOILScores <- allSOILScores[selectedSet]
  selectedSetSorted <- selectedSet[order(selectedSOILScores,decreasing=TRUE)]
  
  numSelected <- length(selectedSet)
  p <- dim(xExp)[2]
  candMat <- matrix(0,nrow=numSelected,ncol=p)
  tempVarSet <- c()
  for (i in 1:numSelected) {
    tempVarSet <- c(tempVarSet,selectedSetSorted[i])
    candMat[i,tempVarSet] <- 1
  }
  
  origCandMat <- candMat
  reRunSOIL_SmallestModel <- SOIL(x=xExp,y=yExp,family="gaussian",weight_type=smallestModelWeightType,
                                  psi=smallestModelPsi,n_train_bound = numModels + 2,
                                  n_train = numModels + 4,
                                  candidate_models = candMat,method="customize")
  
  minInd <- which.max(reRunSOIL_SmallestModel$weight)
  maxInd <- dim(candMat)[1]
  candMat <- candMat[minInd:maxInd,]  
  
  resList <- list(origCandMat = origCandMat,
                  candMat = candMat,
                  minInd = minInd,
                  maxInd = maxInd,
                  reRunSOIL_SmallestModel = reRunSOIL_SmallestModel)
  return(resList)
}