\name{plot.MCRestimate}
\alias{plot.MCRestimate}
\title{Plot method for a objects of class MCRestimate}
\description{plot.MCRestimate visualizes a 'vote matrix'. A 'vote
  matrix' is the result of a classification procedure. For
  every sample (=row) i and every class (=column) j the matrix element
  [i,j] is the probability or frequency the
  classification method predicts sample i as a member of class j.}

\usage{
\method{plot}{MCRestimate}(x,
                 class.factor=NULL,
                 rownames.from.object=FALSE,
                 sample.order=TRUE,
                 legend=FALSE,
                 mypalette=NULL,
                 shading=NULL,
                 xlab="Sample ID",
                 ylab="Frequency of correct classification",
                 cex.axis=1,...)}
\arguments{
  \item{x}{Object of S3 class \code{MCRestimate} or a matrix}
  \item{class.factor}{Factor. Its length must match the number of rows
    in \code{x} and the levels must be the same as the colnames in \code{x}.
    If \code{x} is of class \code{MCRestimate} this argument will be ignored.}
  \item{rownames.from.object}{Logical. If TRUE then the rownames of the
    matrix or the sample names of \code{MCRestimate} in \code{x} are used  as labels for the x-axis}
  \item{sample.order}{Logical. If TRUE then the samples are ordered by
    class membership}
  \item{legend}{Logical. If TRUE then there will be a small legend
    in the output}
  \item{mypalette} {vector with length equal to the number of
    classes. The vector specifies the color for the bar representing the
    classes. If 'NULL' colors choosen by the author are used.}
  \item{shading} {the density of shading lines for the rectagles that indicate the groups, in lines per inch. The default
          value of 'NULL' means that no shading lines are drawn.}
  \item{xlab}{Character}
  \item{ylab}{Character}
  \item{cex.axis}{numeric}
  \item{...}{Further arguments that are passed on to plot.default}
}

\value{The function is called for its side effect, creating a plot on the
  active graphics device.}

\author{Markus Ruschhaupt \url{mailto:m.ruschhaupt@dkfz.de}}

\seealso{\code{\link{MCRestimate}}}


\examples{
  library(MCRestimate)
  x <- c(0.5, 0.3, 0.7, 0.3, 0.8, 0.2, 0.3)
  mat2 <- cbind(x, 1-x)
  fac2 <- factor(c("low", rep("high", 3), rep("low", 3)))
  colnames(mat2) <- levels(fac2)

  mat3 <- cbind(x/3, 2*x/3, 1-x)
  fac3 <- factor(c(rep("high", 3), rep("intermediate", 2), rep("low", 2)))
  colnames(mat3) <- levels(fac3)
if (interactive()) {
  x11(width=9, height=9)
  par(mfrow=c(3,1))}
  plot.MCRestimate(mat2, fac2)
  plot.MCRestimate(mat2, fac2, sample.order=FALSE)
  plot.MCRestimate(mat3, fac3)
}

\keyword{hplot}
