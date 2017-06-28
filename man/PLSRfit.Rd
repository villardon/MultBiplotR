\name{PLSRfit}
\alias{PLSRfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Partial Least Squares Regression (PLSR)
}
\description{
Fits a Partial Least Squares Regression (PLSR) to a set of two continuous data matrices
}
\usage{
PLSRfit(Y, X, S = 2, center = TRUE, scale = TRUE, tolerance = 5e-06, maxiter = 100, show = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
The matrix of dependent variables
}
  \item{X}{
The Matrix of Independent Variables
}
  \item{S}{
Dimension of the solution. The default is 2
}
  \item{center}{
Logical. Should the original data be ccentred.
}
  \item{scale}{
Logical. Should the original data be scaled.
}
  \item{tolerance}{
Tolerance for the algorithm.
}
  \item{maxiter}{
Maximum number of iterations for the algorithm.
}
  \item{show}{
Logical. Should the calculation process be shown on the screen
}
}
\details{
Fits a Partial Least Squares Regression (PLSR) to a set of two continuous data matrices
}
\value{
An object of class "PLSR" 
\item{Method}{PLSR1}
\item{X }{Independent Variables}
\item{Y }{Dependent Variables}
\item{center}{Are data centered?}
\item{scale}{Are data scaled?}
\item{ScaledX}{Scaled Independent Variables}
\item{ScaledY}{Scaled Dependent Variables}
\item{XScores}{Scores for the Independent Variables}
\item{XWeights}{Weights for the Independent Variables - coefficients of the linear combination}
\item{XLoadings}{Factor loadings for the Independent Variables }
\item{YScores}{Scores for the Dependent Variables}
\item{YWeights}{Weights for the Dependent Variables - coefficients of the linear combination}
\item{YLoadings}{Factor loadings for the Dependent Variables }
\item{XStructure}{Structure Correlations for the Independent Variables }
\item{YStructure}{Structure Correlations for the Dependent Variables }
\item{YXStructure}{Structure Correlations two groups}
}
\references{
Wold, S., Sjöström, M., & Eriksson, L. (2001). PLS-regression: a basic tool of chemometrics. Chemometrics and intelligent laboratory systems, 58(2), 109-130.
}
\author{
Jose Luis Vicente Villardon
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
