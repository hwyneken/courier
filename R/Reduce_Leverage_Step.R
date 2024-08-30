Reduce_Leverage_Step <- function(X,y,selectedSet) {
  
  MaxSteps <- 5
  
  currentX <- X
  currentY <- y
  startM <- lm(currentY ~ 0 + currentX[,selectedSet])
  currentHatVals <- hatvalues(startM)
  
  continueMarker <- ifelse(max(currentHatVals) < 0.95,FALSE,TRUE)
  
  currIter <- 0
  while(continueMarker) {
    
    tempX <- currentX[which(currentHatVals < 0.95),]
    tempY <- currentY[which(currentHatVals < 0.95)]
    
    tempM <- lm(tempY ~ 0 + tempX[,selectedSet])
    currIter <- currIter + 1
    continueMarker <- ifelse((max(hatvalues(tempM) < 1) | (currIter >= MaxSteps)),FALSE,TRUE)
    
    currentX <- tempX
    currentY <- tempY
    currentHatVals <- hatvalues(tempM)
  }
  
  if (dim(currentX)[1] <= length(selectedSet)) {
    res <- "Failure"
  }
  else {
    res <- list(X = currentX,y=currentY)
  }
  return(res)
}