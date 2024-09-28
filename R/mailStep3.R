#' mailStep3
#' 
#' \code{mailStep3} Returns a "selected set" of variables based on their SOIL scores \insertCite{ye2018sparsity}{courier}.
#' 
#' @param allSOILScores: The (average over bootstrapped simulation) SOIL scores per variable as returned in step 2. 
#' @param maxModelSize: The maximum model size. When p >= 0.5 * nExp, this is set to min(c(floor(NExp/2), floor(pExp/2))
#' @param xCon: The confirmation data. This is only used to check if the any of the selected variables have a) zero variance in the confirmatory data and b) if there are any linear dependencies in xCon. If a selected variable meets condition a) or b), it will be removed from the selected set.  
#' 
#' @details{
#' This is "Step 3" of the MAIL algorithm.
#' In this step, we create the set of selected variables, often denoted by "A-Hat" in the statistics literature.
#' The most important output of the function is selectedSet.
#' 
#' Let's say that the X matrix is p-dimensional.
#' Therefore, we can represent each variable by its column index.
#' If we were to select all variables, then selectedSet would be 1:p.
#' Note, we only want to return the indices of the selected variables, not their column names.
#' 
#' Unless nExp >= 2*p, we will not select all of the variables. 
#' We will only select variables whose average SOIL score is greater than a cutoff.
#' For more details on how the SOIL scores are calculated, see \code{mailStep2}.
#' The cutoff is determined by the SOIL score with rank maxModelSize.
#' 
#' As an example, consider the case p = 3.
#' Let's say the SOIL scores are 1, 0.2, and 0.5 and maxModelSize = 2.
#' This means that we want to select the two variables with the highest SOIL scores.
#' Thus, the cutoff is 0.5 (2nd highest SOIL score), and selectedSet = c(1,3).
#' 
#' There are two special considerations.
#' First, we could have ties.
#' What if p = 4, allSOILScores = c(1, 0.2, 0.5, 0.5) and maxModelSize = 2.
#' Then, we will make the arbitrary choice to make selectedSet = c(1,3).
#' 
#' Second, we could have selected variables that will mess up linear regression in the confirmatory data.
#' The selected variables could have zero variance in xCon or could be linearly dependent on each other.
#' The function \code{removeUnidentCols} removes these variables from selectedSet.
#' 
#' In principle, selectedSet could be empty, if all variables were unidentifiable in xCon.
#' In that case, selectedSet is integer(0). 
#' The code for MAIL checks to see if length(selectedSet) > 0.
#' If this condition is not met, the algorithm stops.
#' }
#' 
#' @return A list with the following elements:
#' 
#' \itemize{
#'   \item selectedSet: This is the final set of selected variables to be used in the rest of the MAIL algorithm. 
#'   \item origSelectedSet: This is the set of selected variables from the best SOIL scores only - i.e. before checking for dependent columns in xCon. This is for debugging purposes only.
#'   \item numSelected: The number of variables in selectedSet.
#' }
#' 
#' @export

mailStep3 <- function(allSOILScores,maxModelSize,xCon) {
  soilCutoff <- sort(allSOILScores,decreasing=TRUE)[maxModelSize]
  selectedSet <- which(allSOILScores >= soilCutoff)
  
  # for breaking ties
  if (length(selectedSet) > maxModelSize) {
    selectedSet <- selectedSet[order(allSOILScores[selectedSet],decreasing=TRUE)[1:maxModelSize]]
  }
  origSelectedSet <- selectedSet
  
  # for removing variables that are unidentifiable in xCon
  selectedSet <- removeUnidentCols(selectedSet,xCon)
  
  numSelected <- length(selectedSet)
  
  
  resList <- list(selectedSet = selectedSet,
                  origSelectedSet = origSelectedSet,
                  numSelected = numSelected)
}