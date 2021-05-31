\name{CheckBinaryMatrix}
\alias{CheckBinaryMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Checks if a data matrix is binary}
\description{
Checks if a data matrix is binary
}
\usage{
CheckBinaryMatrix(x)
}

\arguments{
  \item{x}{
Matrix to check.
}
}
\details{
Checks if all the entries of the matix are either 0 or 1.
}
\value{
\code{TRUE} if the matrix is binary.
}

\author{
Jose Luis Vicente-Villardon
}

\examples{
\donttest{data(spiders)
sp=Dataframe2BinaryMatrix(spiders)
CheckBinaryMatrix(sp)}

}
