\name{class.factor.format}
\alias{class.factor.format}
\title{A function for creating a factor from the phenoData slot of an exprSet}
\usage{class.factor.format(x, class.column, reference.class=NULL)}

\arguments{
  \item{x}{an \code{exprSet}}
  \item{class.column}{either a number or a  character string specifying
    the relevant column of the phenoData slot in \code{x}}
  \item{reference.class}{character. If specified the result will be
   is a factor with only 2 levels: the reference class versus all other}
}
\description{This function creates a factor whose levels represent the
  different classes for a classification problem. It is derived from the
  column specified in the argument \code{class.column}. If the factor
  has more than two levels, the argument \code{reference.class} can be
  used to transform this into a two-class problem: reference class
  versus the rest.}
\value{A factor of length nrow(pData(x))}
\author{Markus Ruschhaupt \url{mailto:m.ruschhaupt@dkfz.de}}
\examples{library(MCRestimate)
data(exSet)
class.factor.format(exSet, "type")}

\keyword{}
