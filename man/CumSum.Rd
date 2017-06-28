\name{CumSum}
\alias{CumSum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cummulative sums
}
\description{
Cummulative sums 
}
\usage{
CumSum(X, dimens = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Data Matrix
}
  \item{dimens}{
Dimension for summing
}
}
\details{
Cummulative sums within rows (dimens=1) or columns (dimens=2) of a data matrix
}
\value{
A matrix of the same size as X with cummulative sums within each row or each column
}

\author{
Jose Luis Vicente Villardon
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
data(wine)
X=wine[,4:21]
CumSum(X,1)
CumSum(X,2)
}
