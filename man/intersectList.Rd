\name{intersectList}
\alias{intersectList}
\title{A function for creating a all possible intersects for a list of sets.}
\usage{intersectList(x)}

\arguments{
  \item{x}{a \code{list}}
}
\description{The list contains several sets. The function calculates all possible intersections.}
\value{A list containing all possible intersections.}
\author{Markus Ruschhaupt \url{mailto:m.ruschhaupt@dkfz.de}}
\examples{
a <- list (x=1:5,y=3:4,z=c(1,3,3))
intersectList(a)
}

\keyword{file}
