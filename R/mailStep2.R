#' mailStep2
#' 
#' 
#' \code{mailStep2} Runs the SOIL importance calculation for MAIL (2nd step).
#' 
#' @param numSelectionIter: Expects integer >= 1. The number of bootstrapped iterations of SOIL to run. Default = 10 for most algorithms.
#' @param maxModelSize: The maximum model size. 
#' @param xExp: The X matrix of the exploratory set.
#' @param yExp: The y vector of the exploratory set.
#' @param firstSOILWeightType: The weight type choice for SOIL. This can be "AIC", "BIC" or "ARM".
#' @param firstSOILPsi: The complexity penalty for SOIL. This can be any real value - the suggested value is 0.5.
#' @param verbose: If TRUE, print the iteration step (for each iteration in numSelectionIter). (default == FALSE).
#' 
#' @return A list with the following elements:
#' 
#' \itemize{
#'   \item allSOILScores: The average SOIL importance score across numSelectionIter bootstrap iterations.
#'   \item soilRes: If numSelectionIter == 1, then this is just the SOIL output object. If numSelectionIter > 1, this is a list of SOIL outputs for each 
#'        bootstrap iteration. This does not affect the rest of the algorithm, it just allows for transparency and debugging.
#'   \item soilScoreMat: A matrix of dimension numSelectionIter x p, where p is the dimension of xExp. This is full set of SOIL scores for each variable
#'        across all iterations.
#' }
#' 
#' @examples 
#' 
#' This is "Step 2" in the MAIL algorithm, and the first step with any complexity.
#' The purpose of this step is to get SOIL scores for all p variables.
#' The simplest version of this step is to just run SOIL once with the specificed weight type and complexity penalty (psi level).
#' 
#' By default, we chose to use the "union" option of SOIL, which means that for any given dataset, the SOIL scores are determined
#'   by the union of the solution paths of Lasso, MCP and SCAD.
#' Since these variable selection methods can be unstable, we wanted to add a step to add stability to the algorithm.
#' Our strategy is (if the user chooses to do this), to run multiple iterations of SOIL with bootstrapped samples.
#' 
#' If numSelectionIter > 1, then for each iteration, we first get bootstrapped versions of xExp and yExp, and then SOIL is run with these 
#'    data inputs.
#' The SOIL scores for each iteration are kept track of in soilScoreMat.
#' The most important output is allSOILScores, which is the average SOIL score of each variable across all iterations.
#' If numSelectionIter == 1, then we just use the original xExp and yExp.
#' 
#' @references
#' \insertRef{ye2018sparsity}{courier}
#' 
#' @export
#' 
#' 


mailStep2 <- function(numSelectionIter,maxModelSize,xExp,yExp,firstSOILWeightType,firstSOILPsi,verbose=FALSE) {
  p <- dim(xExp)[2]
  NExp <- dim(xExp)[1]
  
  allSOILScores <- rep(0,times=p) ##
  
  ## added bootstrap on 7/29/2024
  if (numSelectionIter > 1) {
    soilScoreMat <- matrix(NA,nrow=numSelectionIter,ncol=p)
    soilRes <- list()
    for (i in 1:numSelectionIter) {
      if (verbose == TRUE) {
        print(sprintf("\tStep 2: Iteration %d",i))
      }
      
      tempInds <- sample(1:NExp,size=NExp,replace=TRUE)
      tempXExp <- xExp[tempInds,]
      tempYExp <- yExp[tempInds]
      if (firstSOILWeightType != "ARM") {
        tempSOILRes <- SOIL(x=tempXExp,y=tempYExp,n_bound = maxModelSize,
                            weight_type=firstSOILWeightType,
                            psi=firstSOILPsi,family="gaussian",method="union")
      }
      else {
        tempSOILRes <- SOIL(x=tempXExp,y=tempYExp,
                            weight_type = "ARM",
                            psi=firstSOILPsi,family="gaussian",method="union",
                            n_train = ceiling(NExp/2)+4)
      }
      soilRes[[i]] <- tempSOILRes
      soilScoreMat[i,] <- as.numeric(tempSOILRes$importance)
      allSOILScores <- allSOILScores + as.numeric(tempSOILRes$importance)
    }
    soilUncertaintyVec <- 1 + apply(soilScoreMat,2,sd,na.rm=TRUE)
    allSOILScores <- allSOILScores / numSelectionIter    
  }
  else {
    if (firstSOILWeightType != "ARM") {
      soilRes <- SOIL(x=xExp,y=yExp,n_bound = maxModelSize,
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
    soilScoreMat <- matrix(soilRes$importance,nrow=1)
    allSOILScores <- allSOILScores + as.numeric(soilRes$importance)    
  }
  
  resList <- list(allSOILScores = allSOILScores,
                  soilRes = soilRes,
                  soilScoreMat = soilScoreMat)
  return(resList)
  
}