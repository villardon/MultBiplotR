\name{CanonicalDistanceAnalysis}
\alias{CanonicalDistanceAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MANOVA and Canonical Analysis of Distances
}
\description{
Performs a MANOVA and a Canonical Analysis based on of Distance Matrices (usally for continuous data)
}
\usage{
CanonicalDistanceAnalysis(Prox, group, dimens = 2, Nsamples = 1000, 
PCoA = "Standard", ProjectInd = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Prox}{
A object containing proximities
}
  \item{group}{
A factor with the group structure of the rows
}
  \item{dimens}{
The dimension of the solution
}
  \item{Nsamples}{
Number of samples for the permutation test. Number of permutations.
}
  \item{PCoA}{
Type of Principal Coordinates for the Canonical Analysis calculated from the Principal coordinates of the Mean Matrix. "Standard" : Standard Principal Coordinates Analysis, "Weighted": Weighted Principal Coordinates Analysis, "WPCA")
}
  \item{ProjectInd}{
Should the individual points be Projected onto the representation For the moment only available for Continuous Data.
}

}
\details{
Performs a MANOVA and a Canonical Analysis based on of Distance Matrices (usally for continuous data). The MANOVA statistics is calculated from a decomposition of the distance matrix based on a factor of a external classification. The significance of the test is calculated using a premutation test. The approach depens only on the distances and then is useful with any kind of data.

The Canonical Representation is calculated from a Principal Coordinates Analysis od the distance matrix among the means. Thus, it is only possible for continuous data. The PCoA representation can be "Standard" using the means directly, "Weighted" weighting each group with its sample size or using weighted Princiopal Components Analysis of the matrix of means.

A measure of the quality of representation of the groups is provided. When possible, the measure is also provided for the individual points.

Soon, a biplot representation will also be developed.
}
\value{
An object of class "CanonicalDistanceAnalysis" with:
\item{Distances }{The Matrix of Distances from which the Analysis has been made}
\item{Groups }{A factor containing the group struture of the individuals}
\item{TSS }{Total sum of squares}
\item{BSS }{Between groups sum of squares}
\item{WSS }{Within groups sum of squares}
\item{Fexp }{Experimental pseudo F-value}
\item{pvalue }{p value based on the permutation test}
\item{Nsamples }{p value based on the permutation test}
\item{ExplainedVariance }{Variances explained by the PCoA}
\item{MeanCoordinates }{Coordinates of the groups for the graphical representation}
\item{Qualities }{Qualities of the representation of the groups}
\item{CummulativeQualities }{Cummulative qualities of the representation of the groups}
\item{RowCoordinates }{Coordinates of the individuals for the graphical representation}
}
\references{
Gower, J. C., & Krzanowski, W. J. (1999). Analysis of distance for structured multivariate data and extensions to multivariate analysis of variance. Journal of the Royal Statistical Society: Series C (Applied Statistics), 48(4), 505-519.

Krzanowski, W. J. (2004). Biplots for multifactorial analysis of distance. Biometrics, 60(2), 517-524.
}
\author{
Jose Luis Vicente Villardon
}
\note{
The MANOVA and the representation of the means can be applied to any Distance althoug the projection of the individuals can be made only for continuous data.
}

\examples{
\donttest{data(iris)
group=iris[,5]
X=as.matrix(iris[1:4])
D=ContinuousProximities(X,  coef = 1)
CDA=CanonicalDistanceAnalysis(D, group, dimens=2)
summary(CDA)
plot(CDA)}
}

