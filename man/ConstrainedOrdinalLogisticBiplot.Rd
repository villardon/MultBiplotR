\name{ConstrainedOrdinalLogisticBiplot}
\alias{ConstrainedOrdinalLogisticBiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constrained Ordinal Logistic Biplot
}
\description{
Constrained Ordinal Logistic Biplot or Redundancy Analysis for Ordinal Data based on logistic responses
}
\usage{
ConstrainedOrdinalLogisticBiplot(Y, X, dim = 2, Scaling = 5, 
tolerance = 1e-05, maxiter = 100, penalization = 0.1, show = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
A binary data matrix
}
  \item{X}{
A matrix of predictors
}
  \item{dim}{
Dimension of the Solution
}
  \item{Scaling}{
Transformation of the columns of the predictor matrix.
}
  \item{tolerance}{
Tolerance for the algorithm
}
  \item{maxiter}{
Maximum number of iterations.
}
  \item{penalization}{
Penalization for the fit (ridge)
}
  \item{show}{
%%     ~~Describe \code{show} here~~
}
}
\details{
Constrained Ordinal Logistic Biplot or Redundancy Analysis for Ordinal Data based on logistic responses.
}
\value{
An ordinal logistic Biplot with the reponse and the predictive variables projected onto it.
}
\references{
Vicente-Villardon, J. L., & Vicente-Gonzalez, L. Redundancy Analysis for Binary Data Based on Logistic Responses in Data Analysis and Rationality in a Complex World. Springer.
}
\author{
Jose Luis Vicente-Villardon
}

\examples{
# not yet
}
