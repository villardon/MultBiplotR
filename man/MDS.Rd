\name{MDS}
\alias{MDS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multidimensional Scaling
}
\description{
Multidimensional Scaling using SMACOF algorithm and Bootstraping the coordinates.
}
\usage{
MDS(Proximities, W = NULL, Model = c("Identity", "Ratio", "Interval", "Ordinal"), 
dimsol = 2, maxiter = 100, maxerror = 1e-06, Bootstrap = FALSE, nB = 200, 
ProcrustesRot = TRUE, BootstrapMethod = c("Sampling", "Permutation"), 
StandardizeDisparities = FALSE, ShowIter = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Proximities}{
An object of class proximities
}
  \item{W}{
A matrix of weigths
}
  \item{Model}{
MDS model. "Identity", "Ratio", "Interval" or "Ordinal".
}
  \item{dimsol}{
Dimension of the solution
}
  \item{maxiter}{
Maximum number of iterations of the algorithm
}
  \item{maxerror}{
Tolerance for convergence of the algorithm
}
  \item{Bootstrap}{
Should Bootstraping be performed?
}
  \item{nB}{
Number of Bootstrap samples.
}
  \item{ProcrustesRot}{
Should the bootstrap replicates be rotated to match the initial configuration using Procrustes?
}
  \item{BootstrapMethod}{
The bootstrap is performed by samplig or permutaing the residuals?
}
  \item{StandardizeDisparities}{
Should the disparities be standardized
}
  \item{ShowIter}{
Show the iteration proccess
}
}
\details{
Multidimensional Scaling using SMACOF algorithm and Bootstraping the coordinates. MDS performs multidimensional scaling of proximity data to find a least- squares representation of the objects in a low-dimensional space.  A majorization algorithm guarantees monotone convergence for optionally transformed, metric and nonmetric data under a variety of models.
}
\value{
An object of class \code{Principal.Coordinates} and  \code{MDS}. The function adds the information of the MDS to the object of class  \code{proximities}. Together with the information about the proximities the object has:
 \item{Analysis}{The type of analysis performed, "MDS" in this case}
 \item{Model}{MDS model used}
 \item{RowCoordinates}{Coordinates for the objects in the MDS procedure}
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
Commandeur, J. J. F. and Heiser, W. J. (1993). Mathematical derivations in the proximity scaling (PROXSCAL) of symmetric data matrices (Tech. Rep. No. RR- 93-03). Leiden, The Netherlands: Department of Data Theory, Leiden University.

Kruskal, J. B. (1964). Nonmetric multidimensional scaling: A numerical method. Psychometrika, 29, 28-42.

De Leeuw, J. & Mair, P. (2009). Multidimensional scaling using majorization: The R package smacof. Journal of Statistical Software, 31(3), 1-30, http://www.jstatsoft.org/v31/i03/

Borg, I., & Groenen, P. J. F. (2005). Modern Multidimensional Scaling (2nd ed.). Springer.

Borg, I., Groenen, P. J. F., & Mair, P. (2013). Applied Multidimensional Scaling. Springer.

Groenen, P. J. F., Heiser, W. J. and Meulman, J. J. (1999). Global optimization in least squares multidimensional scaling by distance smoothing. Journal of Classification, 16, 225-254.

Groenen, P. J. F., van Os, B. and Meulman, J. J. (2000). Optimal scaling by alternating length-constained nonnegative least squares, with application to distance-based analysis. Psychometrika, 65, 511-524.

}

\author{
Jose Luis Vicente Villardon
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{BootstrapSmacof}}
}
\examples{
data(spiders)
Dis=BinaryProximities(spiders)
MDSSol=MDS(Dis, Bootstrap=FALSE)
plot(MDSSol)

}

