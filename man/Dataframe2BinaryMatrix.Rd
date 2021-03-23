\name{Dataframe2BinaryMatrix}
\alias{Dataframe2BinaryMatrix}

\title{Converts a Data Frame into a Binary Data Matrix}
\description{
Converts a Data Frame into a Binary Data Matrix
}
\usage{
Dataframe2BinaryMatrix(dataf, cuttype = "Median", cut = NULL, BinFact = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{
data.frame to be converted
}
  \item{cuttype}{
Type of cut point for continuous variables. Must be "Median" or "Mean".
Does not have any effect for factors
}
  \item{cut}{
Personalized cut value for continuous variables.
}
  \item{BinFact}{
Should I treat a factor with two levels as binary. This means that only a single dummy rather than two is used
}
}
\details{
The function converts a data frame into a Binary Data Matrix  (A matrix with entries either 0 or 1).

Factors with two levels are directly transformed into a column of 0/1 entries.

Factors with more than two levels are converted into a binary submatrix with
as many rows as \code{x} and as many columns as levels or categories. (Indicator matrix)

Integer Variables are treated as factors

Continuous Variables are converted into binary variables using a cut point that can be
the median, the mean or a value provided by the user.
}
\value{
A Binary Data Matrix.
}

\author{
Jose Luis Vicente Villardon
}

\examples{
data(spiders)
Dataframe2BinaryMatrix(spiders)
}
