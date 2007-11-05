\name{PLR}
\alias{PLR}
\alias{predict.PLR}
\title{A function which performs penalised logistic regression
  classification for two groups}
\description{A function which performs penalised logistic regression.}

\usage{PLR(trainmatrix, resultvector, kappa=0, eps=1e-4)
       \method{predict}{PLR}(object,...)}
\arguments{
  \item{resultvector}{a vector which contains the labeling of the samples}
  \item{trainmatrix}{a matrix which includes the data. The rows corresponds to the observations and the colums to the variables.}
  \item{kappa}{value range for penalty parameter. If more that one parameter is specified the one with the lowest AIC will be used.}
  \item{eps} {precision of convergence}
  \item{object}{a fitted PLR model}
  \item{...}{here a data matrix from samples that should be predicted}
}

\details{}

\value{a list with three arguments
  \item{a}{Intercept estimate of the linear predictor}
  \item{b}{vector of estimated regresion coefficients}
  \item{factorlevel}{levels of grouping variable}
  \item{aics}{vector of AIC values with respect to penalty parameter kappa}
  \item{trs}{vector of effective degrees of freedom with respect to penalty parameter kappa}
  }

\author{Axel Benner, Ulrich Mansmann, based on MathLab code by Paul Eilers}

\examples{
library(golubEsets)
data(Golub_Merge)
eSet<-Golub_Merge
X0 <- t(exprs(eSet))
m <- nrow(X0); n <- ncol(X0)
y <- pData(eSet)$ALL.AML
f <- PLR(X0, y,kappa=10^seq(0, 7, 0.5))
if (interactive()) {
  x11(width=9, height=4)
  par(mfrow=c(1,2))
plot(log10(f$kappas), f$aics, type="l",main="Akaike's Information Criterion", xlab="log kappa", ylab="AIC")
plot(log10(f$kappas), f$trs, type="l",xlab="log kappa",
ylab="Dim",main="Effective dimension")
}
}

\keyword{}

