\name{ClassifierBuild}
\alias{ClassifierBuild}
\alias{ClassifierBuild.default}
\title{Building a classifier as a combination of preprocessing and
  classification method}
\description{builds a classifier as a combination of preprocessing and
  classification methods}
\usage{
ClassifierBuild(eset,
                class.column,
                reference.class=NULL,
                classification.fun,
                variableSel.fun ="identity",
                cluster.fun ="identity",
                poss.parameters=list(),
                cross.inner=10,
                rand=123,
                information=TRUE,
                thePreprocessingMethods=c(variableSel.fun,cluster.fun))
}
\arguments{
  \item{eset}{an object of class \code{exprSet} or \code{exprSetRG} }
  \item{class.column}{a number or a character string which indicated the column of the expression set's phenodata containing the class label}
  \item{reference.class}{a character string with the  name of one class - if specified the class will form the first class and all the other classes will form the second class }
 \item{classification.fun}{a character string which names the function that should be used for the classification}
 \item{variableSel.fun}{character string which names the function that should be used for variable selection}
 \item{cluster.fun}{character string which names the function that should be used for clustering the variables}
  \item{thePreprocessingMethods}{vector of character with the names of all preprocessing functions- can be used instead of 'variableSel.fun' and 'cluster.fun' - see details}
 \item{poss.parameters}{a list of possible values for the parameter of the classification method}
 \item{cross.inner}{integer - the number of nearly equal sized parts the train set should be divided into}
 \item{rand}{integer - the random number generator will be put in a reproducible state}
 \item{information}{information - should classificator specific data be given(depends on the wrapper for the classification method)}
}

\value{a \code{list} with the following arguments:
 \item{classifier.for.matrix}{}
 \item{classifier.for.exprSet}{}
 \item{parameter}{a list consisting of the estimated 'best' parameter
   for each cross-validation part}
 \item{class.method}{string which names the function used for the
  classification}
 \item{thePreprocessingMethods}{character string - name of the
   preprocessing functions that have been used}
  \item{cross.inner}{number of blocks for a the inner cross-validation}
  \item{information}{classificator specific data}
}
\author{Markus Ruschhaupt \url{mailto:m.ruschhaupt@dkfz.de}}
\examples{
\dontrun{RF.classifier <- ClassifierBuild(Golub_Train,
+     class.colum, classification.fun = "RF.wrap",
+     thePreprocessingMethods = Preprocessingfunctions,
+     poss.parameters = list.of.poss.parameter,
+     cross.inner = cross.inner)

}
}
\keyword{file}
