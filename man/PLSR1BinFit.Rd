\name{PLSR1BinFit}
\alias{PLSR1BinFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Partial Least Squares Regression with Binary Response
}
\description{
Fits Partial Least Squares Regression with Binary Response
}
\usage{
PLSR1BinFit(Y, X, S = 2, tolerance = 5e-06, maxiter = 100, show =
            FALSE, penalization = 0, cte = TRUE, Algorithm = 2,
            OptimMethod = "CG")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
The response
}
  \item{X}{
The matrix of independent variables
}
  \item{S}{
The Dimension of the solution
}
  \item{tolerance}{
Tolerance for convergence of the algorithm
}
  \item{maxiter}{
Maximum Number of iterations
}
  \item{show}{
Show the steps of the algorithm
}
  \item{penalization}{
Penalization for the Ridge Logistic Regression
}
  \item{cte}{
Should a constant be included in the model?
}
  \item{Algorithm}{
Algorithm used in the calculations
}
  \item{OptimMethod}{
Optimization methods from optimr
}
}
\details{
The procedure uses the algorithm proposed by Bastien et al () to fit a Partial Lest Squares Regression when the response is Binary. The procedure will be later converted into a Biplot to visulize the results.
}
\value{
Still to be finished
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Jose Luis Vicente Villardon
}


\examples{
# No examples yet
}

