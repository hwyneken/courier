LPM_AIC_CV_50Split_Mult = function(XMat,yVec) {
  N = dim(XMat)[1]
  p = dim(XMat)[2]


  NSim = 10 # repeat the whole thing 10 times
  sigma2EstMetaVec = rep(NA,times=NSim)
  for (i in 1:NSim) {
    expInds = sample(1:N,size=floor(N/2),replace=FALSE)

    xExp = XMat[expInds,]
    xCon = XMat[-1*expInds,]

    yExp = yVec[expInds]
    yCon = yVec[-1*expInds]

    NExp = dim(xExp)[1]
    pExp = dim(xExp)[2]

    NCon = dim(xCon)[1]
    pCon = dim(xCon)[2]


    reRunSOIL_AIC = SOIL(x=xExp,y=yExp,family="gaussian",weight_type="AIC",
                         psi=0,method="union",
                         n_bound = floor(NCon*0.5) )
    candMat_Exp = reRunSOIL_AIC$candidate_models_cleaned
    modelWeight_AIC = reRunSOIL_AIC$weight
    largestIndex = which.max(modelWeight_AIC)
    largestModel = which(candMat_Exp[largestIndex,] != 0)

    numCV = 1000
    mseVec = rep(0,times=numCV)
    for (j in 1:numCV) {
      trainInds = sample(1:NCon,size=floor(0.9*NCon),replace=FALSE)
      trainX = xCon[trainInds,largestModel] # which(candMat[maxInd,] != 0)]
      trainY = yCon[trainInds]

      #
      testX = xCon[-1*trainInds,largestModel] #which(candMat[maxInd,] != 0)]
      testY = yCon[-1*trainInds]

      trainDF = data.frame(y = trainY)
      trainDF = cbind(trainDF,trainX)

      tempM = lm(y ~ .,data=trainDF)
      testDF = as.data.frame(testX)
      colnames(testDF) = colnames(trainDF)[-1]
      tempPred = predict(tempM,newdata=testDF)

      mseVec[j] = mean((testY - tempPred)^2,na.rm=TRUE)
    }
    sigma2EstMetaVec[i] = mean(mseVec)
  }
  estSigma2 = mean(sigma2EstMetaVec,na.rm=TRUE)
  return(estSigma2)

}




