\name{EuclideanDistance}
\alias{EuclideanDistance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Classical Euclidean Distance (Pythagorean Distance)
}
\description{
Calculates the eucliden distances among the rows of an euclidean configurations in any dimensions
}
\usage{
EuclideanDistance(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix containing the euclidean configuration
}
}
\details{
eucliden distances among the rows of an euclidean configurations in any dimensions}
\value{
Returns the distance matrix
}

\author{
Jose Luis Vicente Villardon
}

\examples{
x=matrix(runif(20),10,2)
D=EuclideanDistance(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
