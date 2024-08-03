#' mailStep6
#' 
#' @export
#'

mailStep6 <- function(minInd,maxInd,candMat,firstSOILWeightType,firstSOILPsi,xExp,yExp) {
  
  NExp <- dim(xExp)[1]
  
  if (minInd == maxInd) {
    candMat <- matrix(candMat,nrow=1)
    modelWeight <- 1
  }
  else {
    if (firstSOILWeightType != "ARM") {
      finalSOIL <- SOIL(x=xExp,y=yExp,family="gaussian",weight_type=firstSOILWeightType,
                        psi=firstSOILPsi,
                        candidate_models = candMat,method="customize")
    }
    else {
      finalSOIL <- SOIL(x=xExp,y=yExp,family="gaussian",weight_type="ARM",
                        psi=firstSOILPsi,n_train = ceiling(NExp/2)+4,
                        candidate_models = candMat,method="customize")
    }
    modelWeight <- finalSOIL$weight
  }
  
  resList <- list(modelWeight = modelWeight)
  return(resList)
}