#' removUnidentCols
#' 
#' @export


removeUnidentCols <- function(selectedSet,xCon) {
  ## check for dependent columns
  ## start by excluding columns with zero variance in the confirmation set
  sdInCon <- apply(xCon[,selectedSet],2,sd)
  whichHaveZeroSD <- which(sdInCon < 1e-16)
  selectedSet <- setdiff(selectedSet,selectedSet[whichHaveZeroSD]) # remove these
  
  ## now remove columns that are highly unbalanced
  # why do we need this again
  # don't forget to add mcvThresh back in as an argument
  # mostCommonValueFraction <- function(x) {
  #   valCounts <- table(x)
  #   mostCommonInd <- which.max(valCounts)
  #   mostCommonFraction <- valCounts[mostCommonInd] / sum(valCounts,na.rm=TRUE)
  #   return(mostCommonFraction)
  # }
  # mostCommonFractionInCon <- apply(xCon[,selectedSet],2,mostCommonValueFraction)
  # whichAreOverMCVThresh <- which(mostCommonFractionInCon >= mcvThresh)
  # selectedSet <- setdiff(selectedSet,selectedSet[whichAreOverMCVThresh])
  
  # now check for dependent columns
  fullSelectedX <- xCon[,selectedSet] # we just need to make sure that we can fit models on the confirmation set
  fullSelectedRREM <- pracma::rref(fullSelectedX) ## get the reduced row echelon matrix
  # the linearly independent columns start with 1's on the diagonal
  fullSelectedRREM_Diag <- diag(fullSelectedRREM)
  indsToKeep <- which(fullSelectedRREM_Diag == 1)
  selectedSet <- selectedSet[indsToKeep]
  return(selectedSet)
}