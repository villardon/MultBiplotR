\name{Fraction}
\alias{Fraction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Selection of a fraction of the data
}
\description{
Selects a percentage of the data eliminating the observations with higher Mahalanobis distances to the center.
}
\usage{
Fraction(data, confidence = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Two dimensional data set
}
  \item{confidence}{
Percentage to retain. (0-1)
}
}
\details{
The function is used to select a fraction of the data to be plotted for example when clusters are used. The function eliminates the extreme values.
}
\value{
An object of class \code{fraction} with the following fields
\item{data }{The originaldata }
\item{fraction }{The selected data}
\item{confidence }{The percentage selected}

}
\references{
Meulman, J. J., & Heiser, W. J. (1983). The display of bootstrap solutions in multidimensional scaling. Murray Hill, NJ: Bell Laboratories.

Linting, M., Meulman, J. J., Groenen, P. J., & Van der Kooij, A. J. (2007). Stability of nonlinear principal components analysis: An empirical study using the balanced bootstrap. Psychological Methods, 12(3), 359.
}
\author{
Jose Luis Vicente Villardon
}

\seealso{
\code{\link{ConcEllipse}}, \code{\link{AddCluster2Biplot}}
}
\examples{
a=matrix(runif(50), 25,2)
a2=Fraction(a, 0.7)
}

