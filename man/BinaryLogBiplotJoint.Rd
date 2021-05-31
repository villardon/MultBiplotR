\name{BinaryLogBiplotJoint}
\alias{BinaryLogBiplotJoint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Binary logistic biplot with a gradient descent algorithm.
}
\description{
Binary logistic biplot with a gradient descent algorithm.
}
\usage{
BinaryLogBiplotJoint(x, freq = matrix(1, nrow(x), 1), dim = 2, 
ainit = NULL, tolerance = 1e-04, maxiter = 30, penalization = 0.2, 
maxcond = 7, RotVarimax = FALSE, lambda = 0.1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A binary data matrix
}
  \item{freq}{
A vector of frequencies.
}
  \item{dim}{
Dimension of the solution
}
  \item{ainit}{
Initial values for the row coordinates.
}
  \item{tolerance}{
Tolerance
}
  \item{maxiter}{
Maximum number of iterations.
}
  \item{penalization}{
Penalization for the fit (ridge)
}
  \item{maxcond}{
Naximum condition number
}
  \item{RotVarimax}{
Should a Varimax Rotation be used?
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{\dots}{
Aditional arguments
}
}
\details{
Binary logistic biplot with a gradient descent algorithm. Estimates row and column parameters at the same time.
}
\value{
A logistic biplot object.
}
\references{
Vicente-Villardón, J. L., Galindo-Villardón, M. P., & Blázquez-Zaballos, A. (2006). Logistic biplots. Multiple correspondence analysis and related methods. London: Chapman & Hall, 503-521.

Vicente-Villardon, J. L., & Vicente-Gonzalez, L. Redundancy Analysis for Binary Data Based on Logistic Responses in Data Analysis and Rationality in a Complex World. Springer.
}
\author{
Jose Luis Vicente-Villardon
}

\examples{
# not yet
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
