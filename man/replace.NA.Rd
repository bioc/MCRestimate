\name{replace.NA}
\alias{replace.NA}
\title{Replaces in a given numeric matrix NA values per row or per column.}
\description{
Replaces in a given numeric matrix NA values per row or per column.
}
\usage{
replace.NA(x, replacement, byRow = TRUE)
}
\arguments{
  \item{x}{numeric input matrix}
  \item{replacement}{numeric vector containing the values which are used for NA replacement. If byRow = TRUE, this vector must contain as many values as matrix X has rows. Else, this vector must contain as many elements as matrix X has columns.}
  \item{byRow}{logical. If TRUE, then NA values in row n are replaced by the value at position n in the vector
  'replacement'. Else, NA values are replaced according to their column position.}
}

\value{The numeric input matrix with replaced NA values.}
\author{Patrick Warnat \url{mailto:p.warnat@dkfz-heidelberg.de}}

\keyword{file}
