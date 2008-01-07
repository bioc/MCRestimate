\name{MCRindError}
\alias{MCRindError}
\alias{plotIndGroupVotes}
\title{Individual Error of the outer cross-validations}
\description{\code{MCRindError} returns a vector with the individual
  number of incorrect classified samples for each cross-validation}

\usage{
MCRindError(MCRe,perGroup=FALSE)
plotIndGroupVotes(MCRe, PvD= 0.5, dotCol="red", errCol="black", xlab="",ylab="# misclassified samples (mean + SD)",...)
}
                 
\arguments{
  \item{MCRe}{Object of S3 class \code{MCRestimate}}
  \item{perGroup}{logical - should the error per group be calculated as well?}
  \item{PvD}{distance between the error bar and the text}
  \item{dotCol}{color of the dots in the plot}
  \item{errCol}{color of the error bar}
  \item{xlab}{xlab}
  \item{ylab}{ylab}
  \item{...}{further arguments used for plot}
}

\value{ \code{MCRindError} returns a vector of individual errors. If the per group error is also calculated a matrix instead of a vector is returned. plotIndGroupVotes does not return anything.}

\author{Markus Ruschhaupt \url{mailto:m.ruschhaupt@dkfz.de}}

\seealso{\code{\link{MCRestimate}}}


\examples{
library(MCRestimate)
library(golubEsets)
data(Golub_Train)
exSet <- Golub_Train[1:500,]
result1 <-MCRestimate(exSet,"ALL.AML",classification.fun="RF.wrap",cross.outer=3,cross.repeat=3)
MCRindError(result1)
}


\keyword{}


