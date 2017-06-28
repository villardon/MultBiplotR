\name{Three2TwoWay}
\alias{Three2TwoWay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converts a multitable list to a two way matrix 
}
\description{
 Takes a multitable list of matrices X and converts it to a two way matrix with
 the structure required by the Statis programs using a _ to separate variable and occassion
 or study.
}
\usage{
Three2TwoWay(X, whatlines = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The multitable list.
}
  \item{whatlines}{
Concatenate the rows (1) or the columns (2)
}
}
\details{
 Takes a multitable list of matrices X and converts it to a two way matrix with
 the structure required by the Statis programs using a _ to separate variable and occassion
 or study. When whatlines is 1 the final matrix adds the rows of the three dimensional array, then the columns must
 be the same for all studies. When whatlines is 2 the columns are concatenated and then the number
 of rows must be the same for all studies.
}
\value{
A two way matrix
\item{x }{A two way matrix}
}

\author{
Jose Luis Vicente Villardon
}

\examples{
  # No examples yet
}
