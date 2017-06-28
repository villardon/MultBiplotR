\name{diagonal}
\alias{diagonal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diagonal matrix from a vector
}
\description{
Creates a diagonal matrix from a vector
}
\usage{
diagonal(d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
A numerical vector
}
}

\value{
A diagonal matrix wirh the values of vector in the diagonal a zeros elsewhere
}

\author{
Jose Luis Vicente Villardon
}

\examples{
diag(c(1, 2, 3, 4, 5))}
