\name{WeightedPCoA}
\alias{WeightedPCoA}
\title{Weighted Principal Coordinates Analysis}
\usage{
  WeightedPCoA(Proximities, 
  weigths = matrix(1,dim(Proximities$Proximities)[1],1), 
  dimension = 2, tolerance=0.0001)
}
\description{
  Weighted Principal Coordinates Analysis
}
\arguments{
  \item{Proximities}{
A matrix containing the proximities among a set of objetcs
}
  \item{weigths}{
Weigths
}
  \item{dimension}{
Dimension of the solution
}
  \item{tolerance}{
Tolerance for the eigenvalues
}
}

\details{
Weighted Principal Coordinates Analysis
}
\value{data(spiders)
dist=BinaryProximities(spiders)
pco=WeightedPCoA(dist)
An object of class \code{Principal.Coordinates}
}
\references{
Gower, J. C. (2006) Similarity dissimilarity and Distance, measures of. Encyclopedia of Statistical Sciences. 2nd. ed. Volume 12. Wiley

Gower, J.C. (1966). Some distance properties of latent root and vector methods used in multivariate analysis. Biometrika 53: 325-338.

J.R. Demey, J.L. Vicente-Villardon, M.P. Galindo, A.Y. Zambrano, Identifying molecular markers associated with classifications of genotypes by external logistic
biplot, Bioinformatics 24 (2008) 2832.

Cuadras, C. M., Fortiana, J. Metric scaling graphical representation of Categorical Data. Proceedings of Statistics Day, The Center for Multivariate Analysis, Pennsylvania State University, Part 2, pp.1-27, 1995.
}
\author{
Jose Luis Vicente-Villardon
}

\seealso{
\code{\link{BinaryProximities}}
}
\examples{
data(spiders)
dist=BinaryProximities(spiders)
pco=WeightedPCoA(dist)

}


