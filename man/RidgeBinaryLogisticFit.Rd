\name{RidgeBinaryLogisticFit}
\alias{RidgeBinaryLogisticFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fits a binary logistic regression with ridge penalization
}
\description{
This function fits a logistic regression between a dependent variable y and 
  some independent variables x, and solves the separation problem in this type
  of regression using ridge regression and penalization.
}
\usage{
RidgeBinaryLogisticFit(y, xd, freq, tolerance = 1e-05, maxiter = 100, penalization = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
A vector with the values of the dependent variable
}
  \item{xd}{
A matrix with the independent variables
}
  \item{freq}{
Frequencies of each pattern
}
  \item{tolerance}{
Tolerance for the iterations.
}
  \item{maxiter}{
Maximum number of iterations for convergenc~
}
  \item{penalization}{
Penalization used in the diagonal matrix to avoid singularities.
}
}
\details{
Fits a binary logistic regression with ridge penalization
}
\value{
The parameters of the fit
}

\author{
Jose Luis Vicente Villardon
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{RidgeBinaryLogistic}}
}
\examples{
##---- Should be DIRECTLY executable !! ----

}
