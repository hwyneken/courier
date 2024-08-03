#' mailStep3
#' 
#' @export

mailStep3 <- function(allSOILScores,numModels,xCon) {
  soilCutoff <- sort(allSOILScores,decreasing=TRUE)[numModels]
  selectedSet <- which(allSOILScores >= soilCutoff)
  
  if (length(selectedSet) > numModels) {
    selectedSet <- selectedSet[order(allSOILScores[selectedSet],decreasing=TRUE)[1:numModels]]
  }
  origSelectedSet <- selectedSet
  selectedSet <- removeUnidentCols(selectedSet,xCon)
  
  numSelected <- length(selectedSet)
  
  
  resList <- list(selectedSet = selectedSet,
                  origSelectedSet = origSelectedSet,
                  numSelected = numSelected)
}