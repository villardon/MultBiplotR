\name{matrixsqrtinv}
\alias{matrixsqrtinv}

\title{
Inverse of the Matrix squared root
}
\description{
Inverse of the Matrix square root of a matrix using the eigendecomposition.
}
\usage{
matrixsqrtinv(S, tol = sqrt(.Machine$double.eps))
}
\arguments{
  \item{S}{
A squered matrix
}
  \item{tol}{
Tolerance for the igenvalues
}
}
\details{
Inverse of the  Matrix square root of a matrix using the eigendecomposition and removing the eigenvalues smaller than a tolerance
}
\value{
The inverse matrix square root of the argument
}

\author{
Jose Luis Vicente Villardon
}

\seealso{
\code{\link{ginv}}
}
\examples{
data(iris)
x=as.matrix(iris[,1:4])
S= t(x) % * % x
matrixsqrtinv(S)
}
