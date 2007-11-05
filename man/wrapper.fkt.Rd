\name{RF.wrap}
\alias{RF.wrap}
\alias{PAM.wrap}
\alias{PLR.wrap}
\alias{SVM.wrap}
\alias{GPLS.wrap}


\title{Wrapper function for different classification methods}
\description{Wrapper function for different classification methods used
  by \code{MCRestimator}. These functions are mainly used within the
function \code{\link{MCRestimate}}}  

\usage{
RF.wrap(x,y,ntree=500,\dots)
PAM.wrap(x,y,threshold,\dots)
PLR.wrap(x,y,kappa=0,eps=1e-4,...)
SVM.wrap(x,y,gamma = NULL, kernel = "radial", ...)
GPLS.wrap(x,y,\dots)
}

\arguments{
  \item{x,y}{x is a matrix where each row refers to a sample a each
    colum refers to a gene; y is a factor which includes the class for
    each sample}
   \item{threshold}{the threshold for PAM}
   \item{kappa}{the penalty parameter for the penalised logistic regression}
   \item{eps}{}
   \item{ntree}{number of trees to grow in randomForest}
   \item{gamma}{parameter for support vector machines}
   \item{kernel}{parameter for support vector machines}
  \item{\dots}{Further parameters}
}

\value{
  Every function return a predict function which can be used to predict
  the classes for a new data set.}

\author{Markus Ruschhaupt \url{mailto:m.ruschhaupt@dkfz.de}}

\seealso{\code{\link{MCRestimate}}}

\examples{}

\keyword{}
