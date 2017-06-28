\name{GeneralizedProcrustes}
\alias{GeneralizedProcrustes}
\title{Generalized Procrustes Analysis}
\usage{
  GeneralizedProcrustes(x, tolerance = 1e-05, maxiter = 100, Plot = FALSE)
}
\description{
  Generalized Procrustes Analysis
}

\arguments{
  \item{x}{
Three dimensional array with the configurations. The first dimension contains the rows of the configurations,
  the second contains the columns and the third the number of configurations. So x[,,i] is the \emph{i-th}
  configuration
}
  \item{tolerance}{
Tolerance for the Procrustes algorithm. 
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{Plot}{
Should the results be plotted?
}
}

\details{
Generalized Procrustes Analysis for several configurations contained in a three-dimensional
array.
}
\value{
An object of class \code{GenProcustes}.This has components:
\item{History}{History of Iterations}
\item{X}{Initial configurations in a three dimensional array}
\item{RotatedX}{Transformed configurations in a three dimensional array}
\item{Scale}{Scale factors for each configuration}
\item{Rotations}{Rotation Matrices in a three dimensional array}
\item{rss}{Residual Sum of Squares}
\item{Fit}{Goodness of fit as percent of expained variance}
}
\references{
Gower, J.C., (1975). Generalised Procrustes analysis. Psychometrika 40, 33-51.

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
n=dim(spiders)[1]
p=dim(spiders)[2]
prox=array(0,c(n,2,4))

p1=BinaryProximities(spiders,coefficient=5)
prox[,,1]=PrincipalCoordinates(p1)$RowCoordinates
p2=BinaryProximities(spiders,coefficient=2)
prox[,,2]=PrincipalCoordinates(p2)$RowCoordinates
p3=BinaryProximities(spiders,coefficient=3)
prox[,,3]=PrincipalCoordinates(p3)$RowCoordinates
p4=BinaryProximities(spiders,coefficient=4)
prox[,,4]=PrincipalCoordinates(p4)$RowCoordinates
GeneralizedProcrustes(prox)
}
