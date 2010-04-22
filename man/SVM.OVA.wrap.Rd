\name{SVM.OVA.wrap}
\alias{SVM.OVA.wrap}
\title{SVM with 'One-Versus-All' multiclass approach}
\description{
  Multiclass approach where k binary SVM classifiers are constructed for
  a classification problem with k classes:
  Every classifier is trained to distinguish samples of one class from
  samples of all other classes.
  For prediction of the class of a new sample, the sample is classified
  by all k classifiers,
  and the class corresponding to the classifier with the maximum decision value is chosen.
}
\usage{
SVM.OVA.wrap(x,y,gamma = NULL, kernel = "radial", ...)
}
\arguments{
  \item{x,y}{x is a matrix where each row refers to a sample and each column refers to a gene; y is a factor which includes the class for each sample}
  \item{gamma}{parameter for support vector machines}
  \item{kernel}{parameter for support vector machines}
  \item{\dots}{Further parameters}
}
\value{A predict function which can be used to predict the classes for a new data set.}

\author{Patrick Warnat \url{mailto:p.warnat@dkfz-heidelberg.de}}

\seealso{\code{\link{MCRestimate}}}

\examples{
\dontrun{
library(golubEsets)
data(Golub_Train)

class.column <- "ALL.AML"
Preprocessingfunctions <- c("varSel.highest.var")
list.of.poss.parameter <- list(var.numbers = c(250,1000))

Preprocessingfunctions <- c("identity")
class.function <- "SVM.OVA.wrap"
list.of.poss.parameter <- list(gamma = 6)
plot.label <- "Samples"

cross.outer <- 10
cross.repeat <- 20
cross.inner <- 5

SVM.estimate <- MCRestimate(Golub_Train,
		class.column,
		classification.fun = class.function,
		thePreprocessingMethods = Preprocessingfunctions,
		poss.parameters = list.of.poss.parameter,
		cross.outer = cross.outer, cross.inner = cross.inner,
		cross.repeat = cross.repeat, plot.label = plot.label)

}}
\keyword{file}
