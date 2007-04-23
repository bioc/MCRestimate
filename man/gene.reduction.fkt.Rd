\name{varSel.highest.var.eSRG}
\alias{varSel.highest.var.eSRG}
\alias{varSel.green.int.sec.eSRG}
\alias{varSel.green.int.max.eSRG}
\alias{varSel.highest.t.stat.eSRG}
\title{Variable selection and cluster functions to work with MCRestimate.exprSetRG}
\description{Different functions for a variable selection and clustering methods. These functions are used for the function \code{\link{MCRestimate.exprSetRG}}}  

\usage{
       varSel.highest.t.stat.eSRG(sample.gene.matrix,classfactor,theParameter=NULL,var.numbers=500,...)

       varSel.highest.var.eSRG(sample.gene.matrix,classfactor,theParameter=NULL,var.numbers=2000,...)

       varSel.green.int.max.eSRG(sample.gene.matrix,classfactor,theParameter=NULL,lambda=0.5,...)
       varSel.green.int.sec.eSRG(sample.gene.matrix,classfactor,theParameter=NULL, lambda=0.5,...)

}

\arguments{
  \item{sample.gene.matrix}{a matrix in which the rows corresponds to
    genes and the colums corresponds to samples}
  \item{classfactor}{a factor containing the values that should be predicted}
  \item{theParameter}{Either NULL or a logical vector which indicates for every gene if it should
    be left out from further analysis or not} 
  \item{var.numbers}{some methods needs an argument which specifies
      how many variables should be taken}
  \item{lambda}{additional parameter for some methods}
  \item{\dots}{Further parameters}
}

\details{ \code{varSel.highest.var.eSRG} selects a number (specified
  by 'var.numbers') of variables  with the highest variance.
  \code{varSel.highest.t.stat.eSRG} selects variables with highest t
  statistics.  \code{varSel.green.int.max.eSRG} and
  \code{varSel.green.int.sec.eSRG} perform a selection based on the
  green intensity channel. All functions  only work with \code{\link{MCRestimate.exprSetRG}}.} 



  
\value{Every function returns a list consisting of two arguments:
  \item{matrix}{the result matrix of the variable redution or the clustering}
  \item{parameter}{The parameter which are used to reproduce the
    algorithm, i.e. a vector which indicates for every gene if it will
    be left out from further analysis or not if a gene reduction is
    performed or the output of the function kmeans for the clustering
    algorithm.}}

\author{Markus Ruschhaupt \url{mailto:m.ruschhaupt@dkfz.de}}

\seealso{\code{\link{MCRestimate}}}


\examples{library(MCRestimate)
m <- matrix(c(rnorm(10,2,0.5),rnorm(10,4,0.5),rnorm(10,7,0.5),rnorm(10,2,0.5),rnorm(10,4,0.5),rnorm(10,2,0.5)),ncol=2)
cluster.kmeans.mean(m ,number.clusters=3)
}

\keyword{}


