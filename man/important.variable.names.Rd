\name{important.variable.names}
\alias{important.variable.names}
\title{Writing tables with variable information}
\description{The information slot of a
  \code{MCRestimate} object may contain lists
  of variable names that are important for each classification. This
  function produces summary tables of these variables. It is necessary
  that the first column of each information list contains the variable
  names. This is important if someone wants to write a new wrapper for a
  classification method.
}
\usage{
important.variable.names(mcr,file="important_variables",...)}
\arguments{
  \item{mcr}{an object of class \code{MCRestimate} }
  \item{file}{a character string specifies the name of the output files}
   \item{...}{Further arguments that are passed on to plot.default}
}

\value{The function is called for its side effect, writing two tables
  with variable information.}

\author{Markus Ruschhaupt \url{mailto:m.ruschhaupt@dkfz.de}}
\examples{}
\keyword{}
