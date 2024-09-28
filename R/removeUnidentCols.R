#' removUnidentCols
#' 
#' \code{removeUnidentCols} removes variables from the selected set whose full linear targets would not be identifiable in the confirmation data.
#' 
#' @param selectedSet: The set of selected variables which SOIL scores above the cutoff. Selected variables are represented by their index - so if the column names are c("fizz","buzz"), and "fizz" was selected in step 3, then selectedSet = c(1).
#' @param xCon: The X matrix of the confirmatory data. Dimension should be nCon * p.
#' 
#' @details{
#' This function removes two types of unidentifiable variables from selectedSet.
#' The first type is easier to explain - simply put, we check for variables which have practically zero variance
#'  in the confirmation set.
#' In the second step, we check for variables with linearly dependent relationships using the pracma package.
#' 
#' We start the second step by selecting the columns in xCon that correspond to the remaining variables
#'  in selectedSet.
#' Then, we check for linear dependencies.
#' Specifically, we use the \code{pracma::rref} function, which returns a reduced row echelon matrix.
#' 
#' If the data is p-dimensional, the RREF matrix is an upper-triangular matrix. 
#' If there are linear dependencies among the columns, some of the diagonal entries will be zero.
#' A small practical example is given below.
#' In that case, three variables were selected, and the third column depends on the first two.
#' \code{pracma::rref} assigns leading 1's to the first and second columns, but not to the third column.
#' The final step of this function is to filter out all of the variables which do not have a 
#'  diagonal value of "1" in the RREF matrix, and return the filtered selectedSet.  
#' 
#' ```{r}
#' X <- cbind(c(1,1,0),c(0,0,1),c(1,1,1))
#' 
#' X_RREF <- pracma::rref(X)
#' print(X_RREF)
#' ````
#'  
#' }
#' 
#' 
#' @return A (possibly reduced) copy of selected set, with unidentifiable variables removed.
#' 
#' @export


removeUnidentCols <- function(selectedSet,xCon) {
  ## check for dependent columns
  ## start by excluding columns with zero variance in the confirmation set
  sdInCon <- apply(xCon[,selectedSet],2,sd)
  whichHaveZeroSD <- which(sdInCon < 1e-16)
  selectedSet <- setdiff(selectedSet,selectedSet[whichHaveZeroSD]) # remove these
  
  # now check for dependent columns
  fullSelectedX <- xCon[,selectedSet] # we just need to make sure that we can fit models on the confirmation set
  fullSelectedRREM <- pracma::rref(fullSelectedX) ## get the reduced row echelon matrix
  # the linearly independent columns start with 1's on the diagonal
  fullSelectedRREM_Diag <- diag(fullSelectedRREM)
  indsToKeep <- which(fullSelectedRREM_Diag == 1)
  selectedSet <- selectedSet[indsToKeep]
  return(selectedSet)
}