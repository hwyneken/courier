#' MAIL_Split
#'
#' \code{MAIL_Split} runs MAIL with data splitting
#'
#' This is a specific use of the \code{\link{MAIL}} function.
#' The function uses the following arguments with MAIL:
#' \itemize{
#'   \item splitOption = "Split"
#'   \item firstSOILWeightType = "BIC"
#'   \item smallestModelWeightType = "AIC"
#'   \item firstSOILPsi = 0.5
#'   \item smallestModelPsi = 0
#'   sigma2EstFunc = "LPM_AIC_CV_50Split"
#'   verbose = FALSE
#' }
#' @param XMat a n by p numeric matrix
#' @param yVec a n by 1 numeric vector
#' @return resList a list
#'
#' @export
#'
#' @examples
#'
#' #' two main ideas to explain with an example
#' \itemize {
#'   \item how to run MAIL_Split
#'   \item how to run diagnostics for MAIL_Split - is it a reliable method for the current data
#' }
#'
#'
#' #### how to run MAIL_Split
#'
#' # first step - organize the data
#' # gene expression problem: response is the expression of gene probe X249.
#' # install package colonCA through Bioconductor (use these steps)
#' # https://bioconductor.org/packages/release/data/experiment/html/colonCA.html
#' require(colonCA)
#' data(colonCA)
#' dataSet <- t(colonCA@assayData$exprs)
#' colon_y <- dataSet[,249]
#' colon_x <- dataSet[,-1*249]
#' colon_x <- t(unique(t(colon_x))) # remove duplicate columns - original has p = 1999, reduce to p = 1990
#' colon_y <- scale(colon_y)
#' colon_x <- scale(colon_x)
#'
#' colonMAILSplit <- MAIL_Split(XMat = colon_x,yVec = colon_y)
#'
#' @seealso \code{\link{MAIL}} and \code{\link{MAIL_Full}}

MAIL_Split = function(XMat,yVec) {
  resList = MAIL(XMat,yVec,
                 splitOption="Split",
                 firstSOILWeightType = "BIC",
                 smallestModelWeightType = "AIC",
                 firstSOILPsi = 0.5,
                 smallestModelPsi = 0,
                 numSelectionIter = 1,
                 sigma2EstFunc = "LPM_AIC_CV_50Split",
                 trueSD = NULL,
                 verbose=FALSE)
  return(resList)
}
