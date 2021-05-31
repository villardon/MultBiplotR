\name{CheckBinaryVector}
\alias{CheckBinaryVector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Checks if a vector is binary
}
\description{
Checks if all the entries of a vector are 0 or 1
}
\usage{
CheckBinaryVector(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
he vector to check
}
}

\value{
The logical result 
}
\author{
Jose luis Vicente Villardon
}

\examples{
\donttest{x=c(0, 0, 0, 0,  1, 1, 1, 2)
CheckBinaryVector(x)}
}