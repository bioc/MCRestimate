\name{varSel.highest.var}
\alias{varSel.highest.var}
\alias{varSel.highest.t.stat}
\alias{varSel.AUC}
\alias{varSel.removeManyNA}
\alias{varSel.impute.NA}
\alias{identity}
\alias{cluster.kmeans.mean}
\title{Variable selection and cluster functions}
\description{Different functions for a variable selection and clustering methods. These functions are mainly used for the function \code{\link{MCRestimate}}}  

\usage{identity(sample.gene.matrix,classfactor,...)
       varSel.highest.t.stat(sample.gene.matrix,classfactor,theParameter=NULL,var.numbers=500,...)

       varSel.highest.var(sample.gene.matrix,classfactor,theParameter=NULL,var.numbers=2000,...)

       varSel.AUC(sample.gene.matrix, classfactor, theParameter=NULL,var.numbers=200,...)
       cluster.kmeans.mean(sample.gene.matrix,classfactor,theParameter=NULL,number.clusters=500,...)

       varSel.removeManyNA(sample.gene.matrix,classfactor, theParameter=NULL, NAthreshold=0.25,...)
       varSel.impute.NA(sample.gene.matrix ,classfactor,theParameter=NULL,...)
}

\arguments{
  \item{sample.gene.matrix}{a matrix in which the rows corresponds to
    genes and the colums corresponds to samples}
  \item{classfactor}{a factor containing the values that should be predicted}
  \item{theParameter}{Parameter that depends on the function. For
    'cluster.kmeans.mean' either NULL or an output of the function
    \code{kmeans}. If it is NULL then \code{kmeans} will be used to
    form clusters of the genes. Otherwise the already existing clusters
    will be used. In both ways there will be a calculation of the
    metagene intensities afterwards. For the other functions either
    NULL or a logical vector which indicates for every gene if it should
    be left out from further analysis or not}
  \item{number.clusters}{parameter which specifies the number of clusters}  
  \item{var.numbers}{some methods needs an argument which specifies
      how many variables should be taken}
  \item{NAthreshold}{integer- if the percentage of the NA is higher than
    this threshold the variable will be deleted} 
  \item{\dots}{Further parameters}
}

\details{ \code{metagene.kmeans.mean} performs a kmeans clustering with
  a number of clusters specified by 'number clusters' and takes the mean
  of each cluster. \code{varSel.highest.var} selects a number (specified
  by 'var.numbers') of variables  with the highest variance. \code{varSel.AUC} chooses the
  most discriminating variables due to the AUC criterium (the
  library \code{ROC} is required).}



  
\value{Every function returns a list consisting of two arguments:
  \item{matrix}{the result matrix of the variable reduction or the clustering}
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

\keyword{file}


