\name{BinaryLogisticBiplot}
\alias{BinaryLogisticBiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Binary Logistic Biplot
}
\description{
Fits a binary lo gistic biplot to a binary data matrix.
}
\usage{
BinaryLogisticBiplot(x, dim = 2, compress = FALSE, init = "mca", 
method = "EM", rotation = "none", tol = 1e-04, 
maxiter = 100, penalization = 0.2, similarity = "Simple_Matching", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The binary data matrix
}
  \item{dim}{
Dimension of the solution
}
  \item{compress}{
Compress the data before the fitting (not yet implemented)
}
  \item{init}{
Type of initial configuration. ("random", "mirt", "PCoA", "mca")
}
  \item{method}{
 Method to fit the logistic biplot ("EM", "Joint", "mirt", "JointGD", "AlternatedGD", "External", "Recursive")
}
  \item{rotation}{
Rotation of the solution ("none", "oblimin", "quartimin", "oblimax" ,"entropy",  "quartimax", "varimax",  "simplimax" ) see GPARotation
}
  \item{tol}{
Tolerance for the algorithm
}
  \item{maxiter}{
Maximum number of iterations.
}
  \item{penalization}{
Panalization for the different algorithms
}
  \item{similarity}{
Similarity coefficient for the initial configuration or the external model
}
  \item{\dots}{
Any other argument for each particular method.
}
}
\details{
Fits a binary lo gistic biplot to a binary data matrix.

Different Initial configurations can be selected:

1.- random : Random coordinates for each point.

2.- mirt: scores of the procedure mirt (Multidimensional Item Response Theory)

3.- PCoA: Principal Coordinates Analysis

4.- mca: Multiple Correspondence Analysis

We can use also different methods for the estimation

1.- Joint: Joint estimation of the row and column parameters. The Initial alorithm.

2.- EM: Marginal Maximum Likelihood

3.- mirt: Similar to the previous but fitted using the package mirt.

4.- JointGD: Joint estimation of the row and column methods using the gradient descent method.

5.- AlternatedGD: Alternated estimation of the row and column methods using the gradient descent method.

6.- External: Logistic fits on the Principal Coordinates Analysis.

7.- Recursive: Recursive (one axis at a time) estimation of the row and column methods using the gradient descent method. This is similar to the NIPALS algorithm for PCA
}
\value{ 
A Logistic Biplot object.
}
\references{
Vicente-Villardon, J. L., Galindo, M. P. and Blazquez, A. (2006) Logistic Biplots. In Multiple Correspondence Análisis And Related Methods. Grenacre, M & Blasius, J, Eds,  Chapman and Hall, Boca Raton.

Demey, J., Vicente-Villardon, J. L., Galindo, M.P.  AND Zambrano, A. (2008) Identifying Molecular Markers Associated With Classification Of Genotypes Using External Logistic Biplots. Bioinformatics, 24(24): 2832-2838.
}
\author{
Jose Luis Vicente Villardon
}

\seealso{
\code{\link{BinaryLogBiplotJoint}}, \code{\link{BinaryLogBiplotEM}}, \code{\link{BinaryLogBiplotGD}}, \code{\link{BinaryLogBiplotMirt}}, 
}
\examples{
# data(spiders)
# X=Dataframe2BinaryMatrix(spiders)

# logbip=BinaryLogBiplotGD(X,penalization=0.1)
# plot(logbip, Mode="a")
# summary(logbip)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
