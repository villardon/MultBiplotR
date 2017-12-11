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
PLSR1BinFit(Y, X, S = 2, center = TRUE, scale = TRUE, tolerance = 5e-06, maxiter = 100, show = FALSE, penalization = 0, Algorithm=1)
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
  \item{InitTransform}{
Initial Transformation of the data}

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
  \item{Algorithm}{
Algorithm used in the calculations
}
}
\details{
The procedure uses the algorithm proposed by Bastien et al () to fit a Partial Lest Squares Regression when the response is Binary. The procedure will be later converted into a Biplot to visulize the results.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jose Luis Vicente Villardon
}


\examples{
X=as.matrix(wine[,4:21])
colnames(X)= colnames(wine[,4:21])
y=as.numeric(wine[,2])-1
PLSwine=PLSR1BinFit(y, X, S=2,penalization = 0.1)
summary(PLSwine)
}

