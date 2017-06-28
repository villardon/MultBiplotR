\name{ExtractTable}
\alias{ExtractTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extracts unique patterns and its frequencies for a discrete data matrix (numeric)
}
\description{
Extracts the patterns and the frequencies of a discrete data matrix reducing the size of the data matrix in order to accelerate calculations in some techniques.
}
\usage{
ExtractTable(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix of integers containing information of discrete variables. The input matrix must be numerical for the procedure to work properly.
}
}
\details{
For any numerical matrix, calculates the different patterns and the frequencies associated to each pattern
The result contains the pattern matrix, a vector with the frequencies, a list with rows sharing the same pattern.
The final pattern matrix has different ordering than the original matrix.
}
\value{

\item{OriginalNames }{Names before grouping the equal rows}
\item{Patterns}{The reduced table with only unique patterns}
\item{EqualRows}{A list with as many components as unique patterns specifying the original rows with each pattern. That will allow for the reconstruction of the initial matrix}
}

\author{
Jose Luis Vicente-Villardon
}

\examples{
data(spiders)
spidersbin=Dataframe2BinaryMatrix(spiders)
spiderstable=ExtractTable(spidersbin)
}

