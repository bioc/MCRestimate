\name{MCRindError}
\alias{MCRindError}
\alias{plotIndGroupVotes}
\title{Individual Error of the outer cross-validations}
\description{
	\code{MCRindError} returns a vector with the individual number of incorrect classified samples for each cross-validation
  \code{plotIndGroupVotes} plots the individual group votes
}

\usage{
MCRindError(MCRe,
						perGroup=FALSE)
						
plotIndGroupVotes(MCRest,
									PvD= 0.5,
									dotCol="red",
									errCol="black",
									xlab="",
									ylab="# misclassified samples (mean + SD)",
									...)						
}
                 
\arguments{
  \item{MCRe}{Object of S3 class \code{MCRestimate}}
	\item{perGroup}{returns a vector with the individual number of incorrect classified samples for each group}

  \item{MCRest}{Object of S3 class \code{MCRestimate}}
  \item{PvD}{Offset of the text that belongs to a specific point in the plot}
  \item{dotCol}{Color of the dots}
  \item{errCol}{Line color between points}
  \item{xlab}{Label of X-Axis}
  \item{ylab}{Label of Y-Axis}
  \item{...}{Advanced options to the \code{plot} command}

}

 }

\value{ \code{MCRindError} returns a vector of individual errors.}

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


