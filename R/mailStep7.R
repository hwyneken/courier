#' mailStep7
#'
#' @export

mailStep7 <- function(candMat,selectedSet,xCon,yCon,modelWeight,estSigma2) {
  p = dim(xCon)[2]
  numCand = dim(candMat)[1]
  selectedSet = which(candMat[numCand,] != 0)
  numModels = numCand
  numSelected = length(selectedSet)

  tempCoefVec <- rep(0,numSelected)
  tempVarVec <- rep(0,numSelected)

  coefList <- list() # list of coefficients from each submodel
  covMatList <- list() # list of information matrices for each submodel
  for (i in 1:numCand) {
    tempX <- xCon[,which(candMat[i,] != 0)]
    if (sum(candMat[i,] != 0,na.rm=TRUE) == 1) {
      tempX <- matrix(tempX,ncol=1)
    }
    colnames(tempX) = paste("V",which(candMat[i,] != 0),sep="")
    tempDF = data.frame(y=yCon)
    tempDF = cbind(tempDF,tempX)
    tempM = lm(y~.,data=tempDF,tol=1e-16) # set tolerance to all for very dependent columns

    coefList[[i]] <- coef(summary(tempM))
    covMatList[[i]] <- summary(tempM)$cov.unscaled
  }


  for (i in 1:numSelected) {
    tempVar = selectedSet[i]
    tempModelInds = which(candMat[,tempVar] != 0)
    smallestModel = min(tempModelInds)
    tempModelWeight = modelWeight[tempModelInds] / sum(modelWeight[tempModelInds],na.rm=TRUE)
    numTempInds = length(tempModelInds)

    tempCoefVec2 <- rep(0,times=numTempInds)
    tempVarVec2 <- rep(0,times=numTempInds)
    for (j in 1:numTempInds) {
      tempInd = tempModelInds[j]

      if (!(paste0("V",tempVar) %in% rownames(coefList[[tempInd]]))) {
        browser()
      }
      tempCoefVec2[j] <- tempModelWeight[j]*coefList[[tempInd]][paste0("V",tempVar),1]

      tempWeight2 = ifelse(tempInd == numCand,
                           tempModelWeight[j]^2,
                           tempModelWeight[j]^2 + 2*tempModelWeight[j]*sum(tempModelWeight[(j+1):numTempInds],na.rm=TRUE))

      tempCovMat <- covMatList[[tempInd]]
      tempVarVec2[j] <- tempWeight2 * diag(tempCovMat)[paste("V",tempVar,sep="")]
    }

    tempCoefVec[i] <- sum(tempCoefVec2,na.rm=TRUE)
    tempVarVec[i] <- sum(tempVarVec2,na.rm=TRUE)
  }

  tempVarVec <- tempVarVec *  estSigma2

  tempCI <- matrix(0,nrow=numSelected,ncol=2)
  tempCI[,1] <- tempCoefVec - 1.96*sqrt(tempVarVec)
  tempCI[,2] <- tempCoefVec + 1.96*sqrt(tempVarVec)

  betaHatMA <- rep(0,times=p)
  betaHatMA[selectedSet] <- tempCoefVec

  resList = list(betaHatMA = betaHatMA,
                 tempCI = tempCI,
                 margVar = tempVarVec)
}
