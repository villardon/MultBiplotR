\name{ExpandTable}
\alias{ExpandTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Expands a compressed table of patterns and frequencies
}
\description{
Expands a compressed table of patterns and frequencies
}
\usage{
ExpandTable(table)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{
A compressed table of patterns and frequencies
}
}
\details{
To simplify the calculations of some of the algorithms we compress the tables by searching for the distinct patterns and its frequencies. This function recovers the original data. It serves also to assign the corrdinates on the biplot to the original individuals.
}
\value{
A matrix with the original data
}
\author{
Jose Luis Vicente Villardon
}

\examples{
##---- Should be DIRECTLY executable !! ----
}
