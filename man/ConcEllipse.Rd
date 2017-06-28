\name{ConcEllipse}
\alias{ConcEllipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Concentration ellipse for a se of two-dimensional points
}
\description{
The function calculates a non-parametric concentration ellipse for a set of two-dimensional points.
}
\usage{
ConcEllipse(data, confidence=1, npoints=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The set of two-dimensional points
}
  \item{confidence}{
Percentage of points to be included in the ellipse
}
  \item{npoints}{
Number of points to draw the eelipse contour. The hier the number of points the smouther is the ellipse.
}
}
\details{
The procedre uses the Mahalanobis distances to determine the points that will be used for the calculations.
}
\value{
A list with the following fields
\item{data }{Data Used for the calculations}
\item{confidence }{The confidence level used}
\item{ellipse }{The points on the ellipse contour to be plotted}
\item{center }{The center of the points}
}
\references{
Meulman, J. J., & Heiser, W. J. (1983). The display of bootstrap solutions in multidimensional scaling. Murray Hill, NJ: Bell Laboratories.

Linting, M., Meulman, J. J., Groenen, P. J., & Van der Kooij, A. J. (2007). Stability of nonlinear principal components analysis: An empirical study using the balanced bootstrap. Psychological Methods, 12(3), 359.
}
\author{
Jose Luis Vicente Villardon
}

\examples{
data(iris)
dat=as.matrix(iris[1:50,1:2])
plot(iris[,1], iris[,2],col=iris[,5], asp=1)
E=ConcEllipse(dat, 0.95)
plot(E)

}

