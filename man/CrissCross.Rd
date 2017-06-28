\name{CrissCross}
\alias{CrissCross}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Alternated Least Squares Biplot
}
\description{
Alternated Least Squares Biplot with any choice of weigths for each element of the data matrix
}
\usage{
CrissCross(x, w = matrix(1, dim(x)[1], dim(x)[2]), dimens = 2, a0 = NULL, 
b0 = NULL, maxiter = 100, tol = 1e-04, addsvd = TRUE, lambda = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data Matrix to be analysed
}
  \item{w}{
Weights matrix. Must be of the same size as X.
}
  \item{dimens}{
Dimension of the solution.
}
  \item{a0}{
Starting row coordinates. Random coordinates are calculated if the argument is NULL.
}
  \item{b0}{
Starting column coordinates. Random coordinates are calculated if the argument is NULL.
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{tol}{
Tolerance for the algorithm to converge.
}
  \item{addsvd}{
Calculate an additional SVD at the end of the algorithm. That meakes the final solution more readable
}
  \item{lambda}{
Constant to add to the diagonal of the natrices to be inverted in order to improve stability when the matrices are ill-conditioned.
}
}
\details{
The function calculates Alternated Least Squares Biplot with any choice of weigths for each element of the data matrix. The function is useful when we want a low rank approximation of a data matrix in which each element of the matrix has a different weight, for example, all the weights are 1 except for the missing elements that are 0, or to model the logarithms of a frequency table using the frequencies as weights.
}

\value{
An object of class .Biplot" with the following components:
\item{n}{Number of Rows}
\item{p}{Number of Columns}
\item{dim}{Dimension of the Biplot}
\item{EigenValues}{Eigenvalues}
\item{Inertia}{Explained variance (Inertia)}
\item{CumInertia}{Cumulative Explained variance (Inertia)}
\item{RowCoordinates}{Coordinates for the rows}
\item{ColCoordinates}{Coordinates for the columns}
\item{RowContributions}{Contributions for the rows}
\item{ColContributions}{Contributions for the columns}
\item{Scale_Factor}{Scale factor for the traditional plot with points and arrows. The row coordinates are multiplied and the column coordinates divided by that scale factor. The look of the plot is better without changing the inner product. For the HJ-Biplot the scale factor is 1.}

  
}
\references{
GABRIEL, K.R. and ZAMIR, S. (1979). Lower rank approximation of matrices by least squares with any choice of weights. Technometrics, 21: 489-498.
}
\author{
Jose Luis Vicente Villardon
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LogFrequencyBiplot}}
}
\examples{
data(Protein)
X=as.matrix(Protein[,3:11])
X = InitialTransform(X, transform=5)$X
bip=CrissCross(X)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Biplot }
\keyword{ Alternated Least Squares }% __ONLY ONE__ keyword per line
