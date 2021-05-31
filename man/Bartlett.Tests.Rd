\name{Bartlett.Tests}
\alias{Bartlett.Tests}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bartlett tests
}
\description{
Bartlett tests foor the columns of a matrix  and a grouping variable
}
\usage{
Bartlett.Tests(X, groups = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A data frame or a matrix containing several numerical variables
}
  \item{groups}{
A factor with the groups
}
}
\details{
Bartlett tests foor the columns of a matrix  and a grouping variable
}
\value{
A matrix with the tests for each column
}
\references{
Bartlett, M. S. (1937). "Properties of sufficiency and statistical tests". Proceedings of the Royal Statistical Society, Series A 160, 268-282 JSTOR 96803
}
\author{
Jose Luis Vicente Villardon
}

\examples{
\donttest{
data(wine)
Bartlett.Tests(wine[,4:8], groups = wine$Origin)}
}
