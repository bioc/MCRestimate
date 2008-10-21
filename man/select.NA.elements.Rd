\name{select.NA.elements}
\alias{select.NA.elements}
\title{Selects NA values of a given numeric matrix}
\description{
Selects of a given numeric matrix rows or columns containing more NA values than defined by a given threshold.
}
\usage{
select.NA.elements(x, NAthreshold, byRow = TRUE)
}
\arguments{
  \item{x}{numeric input matrix}
  \item{NAthreshold}{numeric value between 0 and 1, determining the allowed percentage of NA values per row or column. Rows or columns containing more NA values are selected.}
  \item{byRow}{logical. If TRUE, then rows with more NA values than determined by the argument threshold are selected, else columns.}
}
\value{
logical vector containing the row or column selection. If argument byRow = TRUE, then value contains as many values as the input matrix contains rows, else it contains as many values as the input matrix contains columns.
}
\author{Patrick Warnat \url{mailto:p.warnat@dkfz-heidelberg.de}}

\keyword{file}
