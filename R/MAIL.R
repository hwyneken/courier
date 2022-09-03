#' MAIL
#'
#' \code{MAIL} runs the Model-Averaged Inferential Learning method under different parameter settings
#'
#' The most important choice is whether or not use data splitting.
#' The advantage of data splitting is to mitigate post selection changes to inference.
#' The advantage of using all of the data is to reduce bias.
#'
#' @param XMat a n by p numeric matrix
#' @param yVec a n by 1 numeric vector
#' @param splitOption Mandatory - can take the values "Full" or "Split"
#' @param firstSOILWeightType Mandatory - can take values "AIC", "BIC" or "ARM"
#' @param smallestModelWeightType Mandatory - can take values "AIC", "BIC" or "ARM"
#' @param firstSOILPsi Mandatory - can take any value in [0,1]
#' @param smallestModelPsi Mandatory - can take any value in [0,1]
#' @param numSelectionIter Optional - defaults to 10, must be an integer >= 1
#' @param sigma2EstFunc Mandatory - this is a string of the function that will estimate the error variance using only XMat and yVec. We recommend using "LPM_AIC_CV_50Split". If the error variance is known, use "trueValue" here.
#' @param trueSD Optional unless "trueValue" has given to the previous argument. This is where the user gives the assumed error standard deviation.
#' @param verbose Optional: default is FALSE - set to TRUE if you want to see printed messages about MAIL's progress.
#'
#' @seealso \code{\link{MAIL_Full}} and \code{\link{MAIL_Split}} for specific versions
#' @seealso \code{\link{LPM_AIC_CV_50Split}} for the recommended variance estimation method
#'
#' @export


MAIL = function(XMat,yVec,
                splitOption,
                firstSOILWeightType,
                smallestModelWeightType,
                firstSOILPsi,
                smallestModelPsi,
                numSelectionIter=10,
                sigma2EstFunc,
                trueSD=NULL,
                verbose=FALSE) {
  # 1) select variables: high vs low value on the right dataset
  # 2) calculate weights
  # 3) construct a candidate matrix

  if (verbose==TRUE) {
    print("Step 1: Organize Data")
  }

  N = dim(XMat)[1]
  p = dim(XMat)[2]

  if (splitOption == "Split") {
    expInds = sample(1:N,size=floor(N/2),replace=FALSE)

    xExp = XMat[expInds,]
    xCon = XMat[-1*expInds,]

    yExp = yVec[expInds]
    yCon = yVec[-1*expInds]

    NExp = dim(xExp)[1]
    pExp = dim(xExp)[2]

    NCon = dim(xCon)[1]
    pCon = dim(xCon)[2]
  }
  else {
    xExp = XMat
    xCon = XMat

    yExp = yVec
    yCon = yVec

    NExp = dim(xExp)[1]
    pExp = dim(xExp)[2]

    NCon = dim(xCon)[1]
    pCon = dim(xCon)[2]
  }

  if (verbose == TRUE) {
    print("Step 2: Run First Model Average")
  }

  allSOILScores = rep(0,times=p) ##
  for (i in 1:numSelectionIter) {
    if (verbose == TRUE) {
      print(sprintf("\tStep 2: Iteration %d",i))
    }

    if (firstSOILWeightType != "ARM") {
      soilRes = SOIL(x=xExp,y=yExp,
                     weight_type=firstSOILWeightType,
                     psi=firstSOILPsi,family="gaussian",method="union")
    }
    else {
      soilRes = SOIL(x=xExp,y=yExp,
                     weight_type = "ARM",
                     psi=firstSOILPsi,family="gaussian",method="union",
                     n_train = ceiling(NExp/2)+4)
    }
    allSOILScores = allSOILScores + as.numeric(soilRes$importance)
  }
  allSOILScores = allSOILScores / numSelectionIter


  if (verbose == TRUE) {
    print("Step 3: Select Variables for the Nested Candidate Set")
  }

  numModels = min(c(floor(dim(xExp)[1]/2),floor(dim(xExp)[2]/2)))
  soilCutoff = sort(allSOILScores,decreasing=TRUE)[numModels]
  selectedSet = which(allSOILScores >= soilCutoff)

  if (length(selectedSet) > numModels) {
    selectedSet = selectedSet[order(allSOILScores[selectedSet],decreasing=TRUE)[1:numModels]]
  }
  origSelectedSet = selectedSet
  selectedSet = removeUnidentCols(selectedSet,xCon)

  numSelected = length(selectedSet)

  if (numSelected > 0) {
    selectedSOILScores = allSOILScores[selectedSet]
    selectedSetSorted = selectedSet[order(selectedSOILScores,decreasing=TRUE)]


    ### create the candidate matrix
    if (verbose == TRUE) {
      print("Step 4: Create Candidate Set, and Smallest Model")
    }

    candMat = matrix(0,nrow=numSelected,ncol=p)
    tempVarSet = c()
    for (i in 1:numSelected) {
      tempVarSet = c(tempVarSet,selectedSetSorted[i])
      candMat[i,tempVarSet] = 1
    }

    origCandMat <- candMat
    reRunSOIL_SmallestModel = SOIL(x=xExp,y=yExp,family="gaussian",weight_type=smallestModelWeightType,
                                   psi=smallestModelPsi,n_train_bound = numModels + 2,
                                   n_train = numModels + 4,
                                   candidate_models = candMat,method="customize")

    minInd = which.max(reRunSOIL_SmallestModel$weight)
    maxInd = dim(candMat)[1]

    if (verbose == TRUE) {
      print("Step 5: Estimate sigma^2")
    }

    ##### Variance Estimation
    if (sigma2EstFunc != "trueValue") {
      tempSigma2Func = get(sigma2EstFunc)
      estSigma2 = tempSigma2Func(XMat,yVec)
    }
    else {
      estSigma2 = trueSD^2
    }

    ### now that we have an estimate of sigma^2,
    ### we can throw out the obviously wrong models that are too large
    ### choose the number of variables as "largestIndex" -
    ### in other words choose min AIC-corrected as the cutoff

    if (verbose == TRUE) {
      print("Step 6: Estimate Final Weights")
    }


    candMat = candMat[minInd:maxInd,]
    if (minInd == maxInd) {
      candMat = matrix(candMat,nrow=1)
      modelWeight = 1
    }
    else {
      if (firstSOILWeightType != "ARM") {
        finalSOIL = SOIL(x=xExp,y=yExp,family="gaussian",weight_type=firstSOILWeightType,
                         psi=firstSOILPsi,
                         candidate_models = candMat,method="customize")
      }
      else {
        finalSOIL = SOIL(x=xExp,y=yExp,family="gaussian",weight_type="ARM",
                         psi=firstSOILPsi,n_train = ceiling(NExp/2)+4,
                         candidate_models = candMat,method="customize")
      }
      modelWeight = finalSOIL$weight
    }



    if (verbose == TRUE) {
      print("Step 7: Get MAIL Estimates and CI's")
    }
    mailOutputs = mailStep7(candMat,selectedSet,xCon,yCon,modelWeight,estSigma2)

  }
  else { # no variables were selected
    selectedSOILScores <- NULL
    selectedSetSorted <- NULL
    soilRes <- NULL

    ### create the candidate matrix
    if (verbose == TRUE) {
      print("Step 4: No Variables Selected, Return Empty Candidate Set")
    }

    candMat <- NULL
    origCandMat <- NULL
    reRunSOIL_SmallestModel <- NULL

    if (verbose == TRUE) {
      print("Step 5: Estimate sigma^2")
    }

    ##### Variance Estimation
    if (sigma2EstFunc != "trueValue") {
      tempSigma2Func = get(sigma2EstFunc)
      estSigma2 = tempSigma2Func(XMat,yVec)
    }
    else {
      estSigma2 = trueSD^2
    }

    if (verbose == TRUE) {
      print("Step 6: No Variables Selected, Do Not Estimate Final Weights")
    }
    modelWeight <- NULL

    if (verbose == TRUE) {
      print("Step 7: No Variables Selected, Do Not Calculate MAIL Estimates and CI's")
    }

    tempVarVec <- NULL
    tempCI <- NULL
    betaHatMA <- rep(0,times=p)
  }




  resList <- list(tempCI = mailOutputs$tempCI,
                  selectedSet = selectedSet,
                  margVar = mailOutputs$margVar,
                  betaHat = mailOutputs$betaHatMA,
                  modelWeight = modelWeight,
                  estSigma2 = estSigma2,
                  candMat = candMat,
                  origSelectedSet = origSelectedSet,
                  reRunSOIL_SmallestModel = reRunSOIL_SmallestModel,
                  origCandMat = origCandMat,
                  origSOILRes = soilRes)

  return(resList)

}
