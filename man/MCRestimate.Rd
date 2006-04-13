\name{MCRestimate}
\alias{MCRestimate}
\alias{MCRestimate.default}
\alias{MCRestimate.exprSetRG}
\alias{print.MCRestimate}
\title{Estimation of misclassification error by cross-validation}
\description{Several repetitions of a cross-validation are performed to
  get 'votes' how stable a method is against different partitions into
  training and test set}
\usage{
MCRestimate(eset,
            class.column,
            reference.class=NULL,
            classification.fun,
            variableSel.fun="identity",
            cluster.fun="identity",
            poss.parameters=list(),
            cross.outer=10,
            cross.repeat=3,
            cross.inner=cross.outer,
            plot.label=NULL,
            rand=123,
            stratify=FALSE,
            information=TRUE,
            thePreprocessingMethods=c(variableSel.fun,cluster.fun))}
\arguments{
  \item{eset}{an object of class \code{exprSet} or class \code{\link[arrayMagic]{exprSetRG-class}}}
  \item{class.column}{a number or a character string which indicated the column of the
expression set's phenodata containing the class label}
  \item{reference.class}{a character string - the name of one class - if
specified, the class will form the first class and all the other classes
will form the second class }
 \item{classification.fun}{character string which names the function that should be used for the classification}
 \item{variableSel.fun}{character string which names the function that should be used for the variable selection}
 \item{cluster.fun}{character string which names the function that
   should be used for the clustering of variables}
 \item{thePreprocessingMethods}{vector of character with the names of all preprocessing functions- can be used instead of 'variableSel.fun' and 'cluster.fun' - see details}
 \item{poss.parameters}{a list of possible values for the parameter of the classification, variable selection, and cluster methods}
 \item{cross.outer}{integer  - the number of nearly equal sized parts the sample set should be divided into (outer cross-validation)}
 \item{cross.repeat}{integer - the number of repetitions of the cross-validation procedure}
  \item{cross.inner}{integer - the number of nearly equal sized parts the train set should be divided into (inner cross-validation)}
  \item{plot.label}{name of one column of the phenodata- if specified, the content of this column will form the labels of the x-axis if the 'votematrix' will be plotted with plot.MCRestimate}
  \item{rand}{integer - the random number generator will be put in a reproducible state}
  \item{stratify}{should a stratified version be used for the cross validation?}
  \item{information}{information - should classificator specific data be given(depends on the wrapper for the classification method)}
}

\value{an object of class \code{MCRestimate} which is a \code{list} with fourteen arguments:
 \item{votes}  {a matrix consisting of the different votes for each sample}
 \item{classes}{the class of each sample}
 \item{table}  {a 'confusion' table, shows the number of 'correct prediction' for each class}
 \item{correct.prediction}{a logical vector - indicates if a sample was
   predicted to be a member of the correct class at least as often as it was predicted to be a member of each other class.}
 \item{correct.class.vote}{vector that contains for every sample the vote for it's correct class}
 \item{parameter}{a list consisting of the estimated 'best' parameter for each cross-validation part}
 \item{class.method}{string which names the function used for the classification}
 \item{thePreprocessingMethods}{character string - name of the preprocessing functions that have been used}
 \item{cross.outer}{number of blocks for a the outer cross-validation}
 \item{cross.repeat}{number of outer cross-validation repetitions}
 \item{cross.inner}{number of blocks for a the inner cross-validation}
 \item{sample.names}{names of the sample}
 \item{information}{classificator specific data (if information is TRUE)}
}
\details{The argument 'thePreprocessingMethods' can be used instead of 'variableSel.fun'
   and 'cluster.fun'. In the first versions of MCRestimate it was only
   possible to have one variable selection and one cluster
   functions. Now it is possible to have more than two functions and the
   ordering is arbitrary, e.g. you can have a variable selection
   function, then a cluster function and then a seccond variable
   selection function.

 If MCRestimate is used with an object of class
 \code{\link[arrayMagic]{exprSetRG-class}}, the preprocessing steps can use the
 green and the red channel seperately but the classification methods
 works with green channel - red channel.

 Note: 'correct prediction' means that a sample was predicted to be
 a member of the correct class at least as often as it was predicted to
 be a member of each other class. So in the two class problem a sample
 is also 'correct' if it has been predicted correctly half of the time. 
}

\author{Markus Ruschhaupt \url{mailto:m.ruschhaupt@dkfz.de}, contributions from Patrick Warnat \url{mailto:p.warnat@dkfz-heidelberg.de}}
\examples{
library(MCRestimate)
library(golubEsets)
data(Golub_Test)
G2 <- Golub_Test[1:500,]
result <- MCRestimate(G2,"ALL.AML",classification.fun="RF.wrap",cross.outer=4,cross.repeat=3)
result
if (interactive()) {
  x11(width=9, height=4)}
plot(result)}
\keyword{}
