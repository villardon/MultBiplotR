\name{ginv}
\alias{ginv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
G inverse
}
\description{
Calculates the g-inverse of a squared matrix using the eigen decomposition and removing the eigenvalues smaller than a tolerance.
}
\usage{
ginv(X, tol = sqrt(.Machine$double.eps))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Matrix to calculate the g-inverse
}
  \item{tol}{
Tolerance. 
}
}
\details{
The function is useful to avoid singularities.
}
\value{
Returns the g-inverse
}

\author{
Jose Luis Vicente Villardon
}

\examples{
data(iris)
x=as.matrix(iris[,1:4])
S= t(x) %*% x 
ginv(S)

}




