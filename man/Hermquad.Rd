\name{Hermquad}
\alias{Hermquad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gauss-Hermite quadrature
}
\description{
Find the Gauss-Hermite abscissae and weights.
}
\usage{
Hermquad(N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
Number of nodes of the quadrature
}
}
\details{
Find the Gauss-Hermite abscissae and weights.
}
\value{

\item{X }{A column vector containing the abscissae.}
\item{W }{A vector containing the corresponding weights.}

}
\references{
Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (1992). Numerical Recipes in C: The Art of Scientific Computing (New York. Cambridge University Press, 636-9.

http://www.mathworks.com/matlabcentral/fileexchange/8836-hermite-quadrature/content/hermquad.m
}
\author{
Jose Luis Vicente Villardon (translated from a Matlab function by  Greg von Winckel) 
)
}

\examples{
Hermquad(10)
}
