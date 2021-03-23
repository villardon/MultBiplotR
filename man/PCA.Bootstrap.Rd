\name{PCA.Bootstrap}
\alias{PCA.Bootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Principal Components Analysis with bootstrap confidence intervals.
}
\description{
Calculates a Principal Components Analysis with bootstrap confidence intervals for its parameters
}
\usage{
PCA.Bootstrap(X, dimens = 2, Scaling = "Standardize columns", B = 1000, type = "np")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The original raw data matrix
}
  \item{dimens}{
Desired dimension of the solution.
}
  \item{Scaling}{
Transformation that should be applied to the raw data.
}
  \item{B}{
Number of Bootstrap samples to draw.
}
  \item{type}{
Type of Bootstrap ("np", "pa", "spper", "spres") 
}
}
\details{
The types of bootstrap used are:
   \itemize{
   \item{"np : "}{Non Parametric}
   \item{"pa : "}{parametric (data is obtained from a Multivariate Normal Distribution)}
   \item{"spper : "}{Semi-parametric Residuals are permutated}
   \item{"spres : "}{Semi-parametric Residuals are resampled}
 }
For the moment, only the non-parametric bootstrap is implemented.

The Principal Components (eigenvectors) are obtained using bootstrap samples.

The Row scotes are obtained projecting the completen data matrix into the bootstrap Principal Components. In this way all the individulas have the same number of replications.
}
\value{

\item{Type }{The type of Bootstrap used}
\item{InitTransform }{Transformation of the raw data}
\item{InitData }{Initial data provided to the function'}
\item{TransformedData }{Transformed Data}
\item{InitialSVD }{Singular value decomposition of the transformed data}
\item{InitScores }{Row Scores for the initial Data}
\item{InitCorr }{Correlation among variables and Principal Components for the Initial Data}
\item{Samples }{Matrix containing the members of the Bootstrap Samples}
\item{EigVal }{Matrix containing the eigenvalues (columns) for each bootstrap sample (columns)}
\item{Inertia }{Matrix containing the proportions of accounted variance (columns) for each bootstrap sample (columns)}
\item{Us }{Three-dimensional array containing the left singular vectors for each bootstrap sample}
\item{Vs }{Three-dimensional array containing the right singular vectors for each bootstrap sample}
\item{As }{Projection of the bootstrap sampled matrix onto the bottstrap principal components}
\item{Bs }{Projection of the bootstrap sampled matrix onto the bottstrap principal coordinates}
\item{Scores }{Projection of the original matrix onto the bootstrap principal components}
\item{Struct }{Correlation of the Initial Variabblñes and the PCs for each bootstrap sample}

}
\references{


Daudin, J. J., Duby, C., & Trecourt, P. (1988). Stability of principal component analysis studied by the bootstrap method. Statistics: A journal of theoretical and applied statistics, 19(2), 241-258.

Chateau, F., & Lebart, L. (1996). Assessing sample variability in the visualization techniques related to principal component analysis: bootstrap and alternative simulation methods. COMPSTAT, Physica-Verlag, 205-210.

Babamoradi, H., van den Berg, F., & Rinnan, Å. (2013). Bootstrap based confidence limits in principal component analysis—A case study. Chemometrics and Intelligent Laboratory Systems, 120, 97-105.

Fisher, A., Caffo, B., Schwartz, B., & Zipunnikov, V. (2016). Fast, exact bootstrap principal component analysis for p> 1 million. Journal of the American Statistical Association, 111(514), 846-860.

}
\author{
Jose Luis Vicente Villardon
}



\seealso{
 \code{\link{PCA.Biplot}}
}
\examples{
\dontrun{X=wine[,4:21]
grupo=wine$Group
rownames(X)=paste(1:45, grupo, sep="-")
pcaboot=PCA.Bootstrap(X, dimens=2, Scaling = "Standardize columns", B=1000)
plot(pcaboot, ColorInd=as.numeric(grupo))
summary(pcaboot)
}
}



