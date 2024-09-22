#' \code{Lasso_CV_Sigma2Est} estimates the error variance in a high dimensional problem
#' using the proposed method from \insertCite{reid2016study}{courier}.
#' The best method in their simulation uses the RSS from the 
#' LASSO coefficients, using 5 or 10-fold CV.
#' 
#' @param XMat a n by p numeric matrix
#' @param yVec a n by 1 numeric vector
#' 
#' @importFrom Rdpack reprompt
#' 
#' @references 
#' \insertRef{reid2016study}{courier}
#' 
#' @export

Lasso_CV_Sigma2Est <- function(XMat,yVec) {
  
  m <- cv.glmnet(x=XMat,y=yVec)
  predVals <- as.numeric(predict(m,newx=XMat))
  
  numSelectedVars = sum(coef(m) != 0) - 1 # exclude intercept
  N <- dim(XMat)[1]
  res <- 1 / (N - numSelectedVars) * sum((yVec - predVals)^2)
  return(res)
}