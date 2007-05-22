\name{MCRconfusion}
\alias{MCRconfusion}
\alias{MCRwrongsamples}
\title{Summary tables for MCRestimate objects}
\description{\code{MCRwrongsamples} returns a matrix with all the samples that have a higher frequency of being predicted as a member of a wrong class than of the correct class for at least one classification method. \code{MCRconfusion} summarizes the result of the vote matices}

\usage{
MCRwrongsamples(x,
         col.names=names(x),
         rownames.from.object=TRUE,
         subgroup=NULL,
         freq=FALSE)

MCRconfusion(x,
             col.names=names(x),
             row.names=NULL)	 }
                 
\arguments{
  \item{x}{List of objects of S3 class \code{MCRestimate}}
  \item{col.names}{Vector of strings used for column names. The length must match the number of objects in \code{x}}
  \item{rownames.from.object}{Logical. If TRUE then the sample names of the
   \code{MCRestimate} object in \code{x} are used as row names}
  \item{subgroup}{Logical. If TRUE then only the samples which belongs
    to the specified group are listed in the table}
  \item{freq}{Logical. If TRUE then the frequency with which each sample in the
    table has been misclassified will be printed.}
  \item{row.names}{Vector of strings used for row names. If not specified the names of the groups are used}
}


\value{ \code{MCRwrongsamples} returns a matrix and \code{MCRconfusion} returns a confusion matrix.}

\author{Markus Ruschhaupt \url{mailto:m.ruschhaupt@dkfz.de}}

\seealso{\code{\link{MCRestimate}}}


\examples{
library(MCRestimate)
library(golubEsets)
data(Golub_Train)
exSet <- Golub_Train[1:500,]
result1 <- MCRestimate(exSet,"ALL.AML",classification.fun="RF.wrap",cross.outer=3,cross.repeat=2)
result2 <- MCRestimate(exSet,"ALL.AML",classification.fun="PAM.wrap",poss.parameters=list(threshold=c(0.5,1)),cross.inner=3,cross.outer=3,cross.repeat=2)
MCRwrongsamples(list(result1,result2),subgroup="AML",col.names=c("Random Forest","PAM"))
MCRconfusion(list(result1,result2),col.names=c("Random Forest","PAM"))
}


\keyword{}


