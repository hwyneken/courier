#' Recommended method for estimating the error variance
#'
#' \code{LPM_AIC_CV_50Split} estimates the error variance with the test set RMSE across 1000 data splits
#'
#' @param XMat a n by p numeric matrix
#' @param yVec a n by 1 numeric vector
#' @param numCV number of cross validation splits (default is 1000)
#'
#' The best way to understand how the function is to break apart the name.
#'
#' \itemize{
#' \item LPM: Largest Plausible Model
#' \item AIC: we choose the largest plausible model by minimizing AIC
#' \item CV: we find the LPM on the training set, and get RMSE from the test set
#' \item 50Split: we use a 50/50 split for the training/test set.
#' }
#'
#' This is very similar to the refitted cross validation method from \insertCite{fan2012variance}{MAIL}
#'
#' @importFrom Rdpack reprompt
#'
#'
#' @references
#' \insertRef{fan2012variance}{MAIL}
#'
#' @export

# help with references here:
# https://cran.r-project.org/web/packages/Rdpack/vignettes/Inserting_bibtex_references.pdf

LPM_AIC_CV_50Split = function(XMat,yVec,
                              numCV = 1000) {
  N = dim(XMat)[1]
  p = dim(XMat)[2]

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


  mseVec = rep(0,times=numCV)
  for (i in 1:numCV) {
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

    mseVec[i] = mean((testY - tempPred)^2,na.rm=TRUE)
  }
  estSigma2 = mean(mseVec)
  return(estSigma2)

}
