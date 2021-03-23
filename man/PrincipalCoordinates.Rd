\name{PrincipalCoordinates}
\alias{PrincipalCoordinates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Principal Coordinates Analysis}
\description{
Principal coordinates Analysis for a matrix of proximities obtained from binary, categorical, continuous or mixed data
}
\usage{
PrincipalCoordinates(Proximities, w = NULL, dimension = 2, tolerance = 1e-04, 
Bootstrap = FALSE, BootstrapType = c("Distances", "Products"), nB = 200, 
ProcrustesRot = TRUE, BootstrapMethod = c("Sampling", "Permutation"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Proximities}{
An object of class \code{proximities}.
}
  \item{w}{
An set of weights.
}
  \item{dimension}{
Dimension of the solution
}
  \item{tolerance}{
Tolerance for the eigenvalues
}
  \item{Bootstrap}{
Should Bootstrap be calculated?
}
  \item{BootstrapType}{
Bootstrap on the residuals of the "distance" or "scalar products" matrix.
}
  \item{nB}{
Number of Bootstrap replications
}

  \item{ProcrustesRot}{
Should each replication be rotated to match the initial solution?
}
  \item{BootstrapMethod}{
The replications are obtained "Sampling" or "Permutating" the residuals.
}
}
\details{
Principal Coordinates Analysis for a proximity matrix previously calculated from a matrix of raw data or directly obsrved proximities.
}
\value{
An object of class \code{Principal.Coordinates}. The function adds the information of the Principal Coordinates to the object of class  \code{proximities}. Together with the information about the proximities the object has:
 \item{Analysis}{The type of analysis performed, "Principal Coordinates" in this case}
 \item{Eigenvalues}{The eigenvalues of the PCoA}
 \item{Inertia}{The Inertia of the PCoA}
 \item{RowCoordinates}{Coordinates for the objects in the PCoA}
 \item{RowQualities}{Qualities of representation for the objects in the PCoA}
  \item{RawStress}{Raw Stress values}
 \item{stress1}{stress formula 1}
 \item{stress2}{stress formula 2}
 \item{sstress1}{sstress formula 1}
 \item{sstress2}{sstress formula 2}
 \item{rsq}{Squared correlation between disparities and distances}
 \item{Spearman}{Spearman correlation between disparities and distances}
 \item{Kendall}{Kendall correlation between disparities and distances}
 \item{BootstrapInfo}{The result of the bootstrap calculations}
}
\references{
Gower, J. C. (2006) Similarity dissimilarity and Distance, measures of. Encyclopedia of Statistical Sciences. 2nd. ed. Volume 12. Wiley

Gower, J.C. (1966). Some distance properties of latent root and vector methods used in multivariate analysis. Biometrika 53: 325-338.

J.R. Demey, J.L. Vicente-Villardon, M.P. Galindo, A.Y. Zambrano, Identifying molecular markers associated with classifications of genotypes by external logistic
biplot, Bioinformatics 24 (2008) 2832.
}
\author{
Jose Luis Vicente-Villardon
}

\seealso{
\code{\link{BinaryProximities}}, \code{\link{BootstrapDistance}}, \code{\link{BootstrapDistance}}, \code{\link{BinaryProximities}}
}
\examples{
data(spiders)
Dis=BinaryProximities(spiders)
pco=PrincipalCoordinates(Dis)
Dis=BinaryProximities(spiders)
pco=PrincipalCoordinates(Dis, Bootstrap=TRUE)

}

