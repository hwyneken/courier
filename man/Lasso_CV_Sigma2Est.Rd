% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lasso_CV_Sigma2Est.R
\name{Lasso_CV_Sigma2Est}
\alias{Lasso_CV_Sigma2Est}
\title{\code{Lasso_CV_Sigma2Est} estimates the error variance in a high dimensional problem
using the proposed method from \insertCite{reid2016study}{courier}.
The best method in their simulation uses the RSS from the
LASSO coefficients, using 5 or 10-fold CV.}
\usage{
Lasso_CV_Sigma2Est(XMat, yVec)
}
\arguments{
\item{XMat}{a n by p numeric matrix}

\item{yVec}{a n by 1 numeric vector}
}
\description{
\code{Lasso_CV_Sigma2Est} estimates the error variance in a high dimensional problem
using the proposed method from \insertCite{reid2016study}{courier}.
The best method in their simulation uses the RSS from the
LASSO coefficients, using 5 or 10-fold CV.
}
\references{
\insertRef{reid2016study}{courier}
}
