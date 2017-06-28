\name{SimpleProcrustes}
\alias{SimpleProcrustes}
\title{Simple Procrustes Analysis}
\usage{
  SimpleProcrustes(X, Y, centre = FALSE)
}

\description{
  Simple Procrustes Analysis for two matrices
}

\arguments{
  \item{X}{
Matrix of the first configuration.
}
  \item{Y}{
Matrix of the second configuration. 
}
  \item{centre}{
Should the matrices be centred before the calculations?
}
}

\details{
Orthogonal Procrustes Analysis for two configurations X and Y. 
The first configuration X is used as a reference and the second, 
Y, is transformed to match the reference as much as possible.
X = s Y T + 1t +E = Z + E

}
\value{
An object of class \code{Procrustes}.This has components:
\item{X}{First Configuration}
\item{Y}{Second Configuration}
\item{Yrot}{Second Configuration after the transformation}
\item{T}{Rotation Matrix}
\item{t}{Translation Vector}
\item{s}{Scale Factor}
\item{rsss}{Residual Sum of Squares}
\item{fit}{Goodness of fit as percent of expained variance}
\item{correlations}{Correlations among the columns of X and Z}
}
\references{
Ingwer Borg, I. & Groenen, P. J.F. (2005). Modern Multidimensional Scaling. Theory and Applications. Second Edition. Springer
}
\author{
Jose Luis Vicente-Villardon
}

\seealso{
\code{\link{PrincipalCoordinates}}
}
\examples{
data(spiders)
}
