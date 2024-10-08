#' MAIL_New
#'
#' \code{MAIL_New} runs the Model-Averaged Inferential Learning method under different parameter settings
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


MAIL_New = function(XMat,yVec,
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
    expInds <- sample(1:N,size=floor(N/2),replace=FALSE)
    
    xExp <- XMat[expInds,]
    xCon <- XMat[-1*expInds,]
    
    yExp <- yVec[expInds]
    yCon <- yVec[-1*expInds]
    
    NExp <- dim(xExp)[1]
    pExp <- dim(xExp)[2]
    
    NCon <- dim(xCon)[1]
    pCon <- dim(xCon)[2]
  }
  else {
    xExp <- XMat
    xCon <- XMat
    
    yExp <- yVec
    yCon <- yVec
    
    NExp <- dim(xExp)[1]
    pExp <- dim(xExp)[2]
    
    NCon <- dim(xCon)[1]
    pCon <- dim(xCon)[2]
  }
  

  
  # We made this change because MAIL was only selecting 20 or so variables in the WEV
  #   example, even when the sample size was very large
  if (pExp >= (1/2)*NExp) {
    maxModelSize <- min(c(floor(NExp/2),
                          floor(pExp/2)))
  }
  else {
    maxModelSize <- pExp
  }
  
  if (verbose == TRUE) {
    print("Step 2: Run First Model Average")
  }
  
  mailStep2Res <- mailStep2(numSelectionIter,maxModelSize,xExp,yExp,firstSOILWeightType,firstSOILPsi,verbose=FALSE)
  allSOILScores <- mailStep2Res$allSOILScores

  
  
  
  if (verbose == TRUE) {
    print("Step 3: Select Variables for the Nested Candidate Set")
  }
  mailStep3Res <- mailStep3(allSOILScores,maxModelSize,xCon)
  selectedSet <- mailStep3Res$selectedSet
  numSelected <- mailStep3Res$numSelected
  
  if (numSelected > 0) {
    
    ### filter out rows of data with very high leverage - test change 8/28/2024
    # fullX_Exp_Selected <- xExp[,selectedSet]
    # fullX_Con_Selected <- xCon[,selectedSet]
    # 
    # fullM_Exp <- lm(yExp ~ 0 + fullX_Exp_Selected)
    # fullM_Con <- lm(yCon ~ 0 + fullX_Con_Selected)
    # 
    # veryHighLeverage_Exp <- which(hatvalues(fullM_Exp) >= 0.95)
    # veryHighLeverage_Con <- which(hatvalues(fullM_Con) >= 0.95)
    # 
    # if (length(veryHighLeverage_Exp) > 0) {
    #   xExp <- xExp[-1 * veryHighLeverage_Exp,]
    #   yExp <- yExp[-1 * veryHighLeverage_Exp]
    # }
    # 
    # if (length(veryHighLeverage_Con) > 0) {
    #   xCon <- xCon[-1 * veryHighLeverage_Con,]
    #   yCon <- yCon[-1 * veryHighLeverage_Con]
    # }
    
    ## add function to do this 8/29/2024
    # absCorFunc <- function(x) {abs(cor(x,yExp))}
    # absCorValues <- apply(xExp,2,absCorFunc)
    # absCor_CutoffPoint <- sort(absCorValues)[floor(NExp/2)]
    
    cleanedData_Exp <- Reduce_Leverage_Step(xExp,yExp,selectedSet)
    cleanedData_Con <- Reduce_Leverage_Step(xCon,yCon,selectedSet)
    
    if ((class(cleanedData_Exp) == "character") | (class(cleanedData_Con) == "character")) {
      stop("Data Leverage Cleaning Failed")
      
    }
    else {
      xExp <- cleanedData_Exp$X
      yExp <- cleanedData_Exp$y
      
      xCon <- cleanedData_Con$X
      yCon <- cleanedData_Con$y
    }
    
    ### create the candidate matrix
    if (verbose == TRUE) {
      print("Step 4: Create Candidate Set, and Smallest Model")
    }
    
    mailStep4Res <- mailStep4(allSOILScores,selectedSet,xExp,yExp,smallestModelWeightType,smallestModelPsi,maxModelSize)
    candMat <- mailStep4Res$candMat
    origCandMat <- mailStep4Res$origCandMat
    minInd <- mailStep4Res$minInd
    maxInd <- mailStep4Res$maxInd
    reRunSOIL_SmallestModel <- mailStep4Res$reRunSOIL_SmallestModel
    selectedSetSorted <- mailStep4Res$selectedSetSorted
    
    if (verbose == TRUE) {
      print("Step 5: Estimate sigma^2")
    }
    
    ##### Variance Estimation
    if (sigma2EstFunc != "trueValue") {
      tempSigma2Func <- get(sigma2EstFunc)
      estSigma2 <- tempSigma2Func(XMat,yVec)
    }
    else {
      estSigma2 <- trueSD^2
    }
    
    ### now that we have an estimate of sigma^2,
    ### we can throw out the obviously wrong models that are too large
    ### choose the number of variables as "largestIndex" -
    ### in other words choose min AIC-corrected as the cutoff
    
    if (verbose == TRUE) {
      print("Step 6: Estimate Final Weights")
    }
    mailStep6Res <- mailStep6(minInd,maxInd,candMat,smallestModelWeightType,smallestModelPsi,xExp,yExp)
    modelWeight <- mailStep6Res$modelWeight
    if (any(is.nan(modelWeight) | is.na(modelWeight))) {
      stop("Failed estimation of model weights")
    }
    
    
    
    if (verbose == TRUE) {
      print("Step 7: Get MAIL Estimates and CI's")
    }
    mailStep7_SafeFunc <- safely(function() {mailStep7(candMat,selectedSet,xCon,yCon,modelWeight,estSigma2)})
    mailOutputs <- mailStep7_SafeFunc()$result
    if (is.null(mailOutputs)) {
      stop("Step 7 Failed - Investigate Manually")
    }
    
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
  
  
  dataList <- list(xExp = xExp,xCon = xCon,
                   yExp = yExp,yCon = yCon)
  
  resList <- list(CIMatrix = mailOutputs$tempCI,
                  selectedSet = selectedSet,
                  margVar = mailOutputs$margVar,
                  betaHat = mailOutputs$betaHatMA,
                  modelWeight = modelWeight,
                  estSigma2 = estSigma2,
                  candMat = candMat,
                  origSelectedSet = mailStep3Res$origSelectedSet,
                  reRunSOIL_SmallestModel = reRunSOIL_SmallestModel,
                  origCandMat = origCandMat,
                  origSOILRes = mailStep2Res$soilRes,
                  dataList = dataList,
                  selectedSetSorted = selectedSetSorted,
                  allSOILScores = allSOILScores,
                  soilScoreMat = mailStep2Res$soilScoreMat)
  
  return(resList)
  
}
