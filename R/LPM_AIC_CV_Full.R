LPM_AIC_CV_Full = function(XMat,yVec) {
  N = dim(XMat)[1]
  p = dim(XMat)[2]
  
  reRunSOIL_AIC = SOIL(x=XMat,y=yVec,family="gaussian",weight_type="AIC",
                       psi=0,method="union",
                       n_bound = floor(N*0.5))
  candMat_Exp = reRunSOIL_AIC$candidate_models_cleaned
  modelWeight_AIC = reRunSOIL_AIC$weight
  largestIndex = which.max(modelWeight_AIC)
  largestModel = which(candMat_Exp[largestIndex,] != 0)
  
  numCV = 1000
  mseVec = rep(0,times=numCV)
  for (i in 1:numCV) {
    trainInds = sample(1:N,size=floor(0.9*N),replace=FALSE)
    trainX = XMat[trainInds,largestModel] # which(candMat[maxInd,] != 0)]
    trainY = yVec[trainInds]
    
    #
    testX = XMat[-1*trainInds,largestModel] #which(candMat[maxInd,] != 0)]
    testY = yVec[-1*trainInds]
    
    trainDF = data.frame(y = trainY)
    trainDF = cbind(trainDF,trainX)
    
    tempM = lm(y ~ .,data=trainDF)
    testDF = as.data.frame(testX)
    colnames(testDF) = colnames(trainDF)[-1]
    tempPred = predict(tempM,newdata=testDF)
    
    mseVec[i] = mean((testY - tempPred)^2,na.rm=TRUE)
  }
  estSigma2 = mean(mseVec)
  return(estSigma2)
}