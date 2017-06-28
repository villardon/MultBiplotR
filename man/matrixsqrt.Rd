\name{matrixsqrt}
\alias{matrixsqrt}

\title{
Matrix squared root
}
\description{
Matrix square root of a matrix using the eigendecomposition.
}
\usage{
matrixsqrt(S, tol = sqrt(.Machine$double.eps))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
A squered matrix
}
  \item{tol}{
Tolerance for the igenvalues
}
}
\details{
Matrix square root of a matrix using the eigendecomposition and removing the eigenvalues smaller than a tolerance
}
\value{
The matrix square root of the argument
}

\author{
Jose Luis Vicente Villardon
}

\examples{
data(iris)
x=as.matrix(iris[,1:4])
S= t(x) % * % x
matrixsqrt(S)
}

